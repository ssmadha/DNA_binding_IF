import gzip
import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmcrameri.cm as cmc
from matplotlib.cm import ScalarMappable
from matplotlib.colors import ListedColormap
from sklearn.cluster import KMeans

# ============================================================
# Settings
# ============================================================

# All genes, k-means clustered (genes and tissues), with k chosen by the
# elbow method over these ranges
GENE_K_VALUES = range(3, 50)
TISSUE_K_VALUES = range(2, 21)
KMEANS_SEED = 0
# Leave the lowest-ASIF gene cluster (mostly ASIF ~ 0) out of the heatmap;
# clustering and the cluster tables still use all genes
DROP_LOWEST_GENE_CLUSTER = False
# Figure size in inches; rows are squeezed to fit the height
FIG_SIZE = (12, 14)
# Minimum PNG resolution; raised automatically so every gene row gets at
# least one pixel
HEATMAP_DPI = 600
FONT_SIZE = 15
plt.rcParams["font.size"] = FONT_SIZE
# Shorter display names for some tissues
TISSUE_LABELS = {"parathyroid gland": "parathyroid"}

# One heatmap per gene cluster, with gene names. The lowest-ASIF cluster
# (7,639 genes, mostly ASIF ~ 0) is left out.
PLOT_LOWEST_CLUSTER_ALONE = False
CLUSTER_ROW_INCHES = 0.15
CLUSTER_GENE_FONT_SIZE = 10
CLUSTER_FIG_WIDTH = 12
CLUSTER_DPI = 300
# Acrobat can't open PDFs taller than 200 in, and Agg can't draw images
# larger than 2^16 px on a side
MAX_PDF_INCHES = 200
MAX_PNG_PIXELS = 65000


# ============================================================
# File paths
# ============================================================

ppi_file = "results_segment_ppi/expressed_coding_isoforms_with_relative_tpm_threshold_1_ASIF.tsv"
# Gene names come from the same Ensembl GTF the pipeline uses
gtf_file = "reference_data/Homo_sapiens.GRCh38.109.gtf.gz"

output_dir = "ASIF_PPI_Figures"
os.makedirs(output_dir, exist_ok=True)


# ============================================================
# 1. Read files
# ============================================================

ppi_df = pd.read_csv(ppi_file, sep="\t")


# ============================================================
# 2. Gene-level ASIF
# ============================================================

asif_columns = [
    col for col in ppi_df.columns
    if col.endswith("_asif")
]

# Take maximum ASIF across transcripts for each gene
all_gene_asif = ppi_df.groupby(
    "gene_id"
)[asif_columns].max()

all_gene_asif.columns = [
    TISSUE_LABELS.get(col.replace("_asif", ""), col.replace("_asif", ""))
    for col in all_gene_asif.columns
]


# ============================================================
# 3. Gene names
# ============================================================

gene_names = {}
with gzip.open(gtf_file, "rt") as handle:
    for line in handle:
        fields = line.split("\t")
        if len(fields) > 8 and fields[2] == "gene":
            gene_id = re.search(r'gene_id "([^"]+)"', fields[8]).group(1)
            name = re.search(r'gene_name "([^"]+)"', fields[8])
            if name:
                gene_names[gene_id] = name.group(1)


# ============================================================
# 4. k-means clustering (elbow method)
# ============================================================

def elbow_kmeans(values, k_values, name):
    """Fit k-means for each k, pick the elbow, and save the elbow plot.

    The elbow is the k whose (k, inertia) point lies farthest from the
    straight line joining the first and last points of the curve, with
    both axes scaled to 0-1. Returns the k-means fit at the elbow.
    """
    fits = {
        k: KMeans(n_clusters=k, n_init=10, random_state=KMEANS_SEED).fit(values)
        for k in k_values
    }
    ks = np.array(list(fits))
    inertias = np.array([fits[k].inertia_ for k in ks])

    x = (ks - ks[0]) / (ks[-1] - ks[0])
    y = (inertias - inertias[-1]) / (inertias[0] - inertias[-1])
    # Distance from the line y = 1 - x
    distances = np.abs(x + y - 1) / np.sqrt(2)
    best_k = int(ks[np.argmax(distances)])

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(ks, inertias, marker="o", linestyle="-")
    ax.axvline(best_k, color="grey", linestyle="--")
    ax.annotate(
        f"k = {best_k}",
        xy=(best_k, fits[best_k].inertia_),
        xytext=(8, 8),
        textcoords="offset points"
    )
    ax.set_title(f"Elbow Plot ({name})")
    ax.set_xlabel("Number of Clusters (k)")
    ax.set_ylabel("Distortion (Inertia)")
    ax.grid(True)
    fig.savefig(
        os.path.join(output_dir, f"PPI_ASIF_kmeans_elbow_{name}.png"),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close(fig)

    print(f"Elbow for {name}: k = {best_k}")
    return fits[best_k]


def cluster_order(values, labels):
    """Order clusters by mean ASIF (highest first), members likewise."""
    values = pd.DataFrame(values)
    row_mean = values.mean(axis=1).to_numpy()
    cluster_means = pd.Series(row_mean).groupby(labels).mean()
    rank = cluster_means.rank(ascending=False, method="first").astype(int) - 1
    cluster_rank = rank.loc[labels].to_numpy()
    order = np.lexsort((-row_mean, cluster_rank))
    return order, cluster_rank


gene_kmeans = elbow_kmeans(all_gene_asif.to_numpy(), GENE_K_VALUES, "genes")
tissue_kmeans = elbow_kmeans(all_gene_asif.T.to_numpy(), TISSUE_K_VALUES, "tissues")

gene_order, gene_clusters = cluster_order(
    all_gene_asif.to_numpy(), gene_kmeans.labels_
)
tissue_order, tissue_clusters = cluster_order(
    all_gene_asif.T.to_numpy(), tissue_kmeans.labels_
)

clustered_values = all_gene_asif.iloc[gene_order, tissue_order]

# Cluster assignments, in figure order (clusters numbered from the top)
pd.DataFrame({
    "gene_id": clustered_values.index,
    "gene_name": [gene_names.get(g, g) for g in clustered_values.index],
    "cluster": gene_clusters[gene_order] + 1,
    "mean_asif": clustered_values.mean(axis=1).round(4).values,
}).to_csv(
    os.path.join(output_dir, "PPI_ASIF_kmeans_gene_clusters.tsv"),
    sep="\t",
    index=False
)
pd.DataFrame({
    "tissue": clustered_values.columns,
    "cluster": tissue_clusters[tissue_order] + 1,
}).to_csv(
    os.path.join(output_dir, "PPI_ASIF_kmeans_tissue_clusters.tsv"),
    sep="\t",
    index=False
)


# ============================================================
# 5. Clustered heatmap
# ============================================================

def cluster_cmap(n_clusters):
    """Categorical colours for cluster strips (tab20, then tab20b)."""
    colors = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
    return ListedColormap([colors[i % len(colors)] for i in range(n_clusters)])


n_gene_clusters = gene_kmeans.n_clusters
n_tissue_clusters = tissue_kmeans.n_clusters

ordered_gene_clusters = gene_clusters[gene_order]
shown_genes = np.ones(len(gene_order), dtype=bool)
if DROP_LOWEST_GENE_CLUSTER:
    # Clusters are ranked by mean ASIF, so the last one is the lowest
    shown_genes = ordered_gene_clusters != n_gene_clusters - 1

shown_values = clustered_values.to_numpy()[shown_genes]
shown_gene_clusters = ordered_gene_clusters[shown_genes]
n_rows, n_cols = shown_values.shape

print(
    f"k-means heatmap: {n_rows:,} of {len(gene_order):,} genes"
    + (" (lowest cluster left out)" if DROP_LOWEST_GENE_CLUSTER else "")
)

fig = plt.figure(figsize=FIG_SIZE)

# Cluster strips on the left (genes) and top (tissues) of the heatmap
gs = fig.add_gridspec(
    nrows=2,
    ncols=2,
    width_ratios=[1, 30],
    height_ratios=[1, 30],
    wspace=0.02,
    hspace=0.02
)

ax_tissue_strip = fig.add_subplot(gs[0, 1])
ax_gene_strip = fig.add_subplot(gs[1, 0])
ax_heatmap = fig.add_subplot(gs[1, 1])

heatmap_height_inches = ax_heatmap.get_position().height * FIG_SIZE[1]
png_dpi = max(HEATMAP_DPI, int(np.ceil(n_rows / heatmap_height_inches)))

image = ax_heatmap.imshow(
    shown_values,
    aspect="auto",
    interpolation="nearest",
    cmap=cmc.batlow,
    vmin=0,
    vmax=1
)
ax_heatmap.set_xticks(range(n_cols))
ax_heatmap.set_xticklabels(clustered_values.columns, rotation=90)
ax_heatmap.set_yticks([])
ax_heatmap.set_xlabel("Tissue")

# Gene cluster strip on the left, with cluster numbers
ax_gene_strip.imshow(
    shown_gene_clusters[:, None],
    aspect="auto",
    interpolation="nearest",
    cmap=cluster_cmap(n_gene_clusters),
    vmin=-0.5,
    vmax=n_gene_clusters - 0.5
)
ax_gene_strip.set_xticks([])
ax_gene_strip.set_yticks([])
boundaries = np.flatnonzero(np.diff(shown_gene_clusters)) + 1
starts = np.concatenate(([0], boundaries))
ends = np.concatenate((boundaries, [n_rows]))
centres = (starts + ends - 1) / 2

# Cluster numbers next to the strip, pushed apart (with a leader line to
# the cluster) where small clusters would make them overlap
min_gap_rows = 1.2 * FONT_SIZE / 72 / heatmap_height_inches * n_rows
label_rows = centres.copy()
for i in range(1, len(label_rows)):
    label_rows[i] = max(label_rows[i], label_rows[i - 1] + min_gap_rows)

strip_transform = ax_gene_strip.get_yaxis_transform()
for cluster, centre, label_row in zip(
    shown_gene_clusters[starts], centres, label_rows
):
    ax_gene_strip.annotate(
        str(cluster + 1),
        xy=(0, centre),
        xycoords=strip_transform,
        xytext=(-1.2, label_row),
        textcoords=strip_transform,
        ha="right",
        va="center",
        arrowprops=dict(arrowstyle="-", linewidth=0.8, shrinkA=2, shrinkB=0)
    )
ax_gene_strip.set_ylabel(
    f"Genes (n = {n_rows:,}), k-means clusters",
    labelpad=3 * FONT_SIZE
)

# Tissue cluster strip on top, with the tissue names above it
ax_tissue_strip.imshow(
    tissue_clusters[tissue_order][None, :],
    aspect="auto",
    interpolation="nearest",
    cmap=cluster_cmap(n_tissue_clusters),
    vmin=-0.5,
    vmax=n_tissue_clusters - 0.5
)
ax_tissue_strip.set_yticks([])
ax_tissue_strip.set_xticks(range(n_cols))
ax_tissue_strip.set_xticklabels(clustered_values.columns, rotation=90)
ax_tissue_strip.xaxis.set_ticks_position("top")
ax_tissue_strip.xaxis.set_label_position("top")

# The axes title sits above the tissue names
ax_tissue_strip.set_title(
    f"ASIF, k-means Clustered\n"
    f"({n_gene_clusters} gene / {n_tissue_clusters} tissue clusters)",
    fontsize=FONT_SIZE + 2
)

for extension in ["png", "pdf"]:
    fig.savefig(
        os.path.join(output_dir, f"PPI_ASIF_all_genes_kmeans_heatmap.{extension}"),
        dpi=png_dpi,
        bbox_inches="tight"
    )
plt.close(fig)


# ============================================================
# 6. ASIF colour legend, saved on its own
# ============================================================

fig, ax_cbar = plt.subplots(figsize=(4, 0.4))
colorbar = fig.colorbar(
    ScalarMappable(norm=image.norm, cmap=image.cmap),
    cax=ax_cbar,
    orientation="horizontal",
    ticks=np.linspace(0, 1, 6)
)
fig.savefig(
    os.path.join(output_dir, "PPI_ASIF_all_genes_kmeans_heatmap_colorbar.png"),
    dpi=300,
    bbox_inches="tight"
)
plt.close(fig)


# ============================================================
# 7. Each gene cluster on its own, with gene names
# ============================================================

def gene_labels(gene_ids):
    """Gene names, keeping the Ensembl ID for missing or repeated names."""
    labels = pd.Series(
        [gene_names.get(gene_id, gene_id) for gene_id in gene_ids],
        index=gene_ids
    )
    duplicated = labels.duplicated(keep=False)
    labels[duplicated] = [
        f"{label} ({gene_id})"
        for gene_id, label in labels[duplicated].items()
    ]
    return labels.values


cluster_dir = os.path.join(output_dir, "kmeans_clusters")
os.makedirs(cluster_dir, exist_ok=True)

clusters_to_plot = range(n_gene_clusters)
if not PLOT_LOWEST_CLUSTER_ALONE:
    # Clusters are ranked by mean ASIF, so the last one is the lowest
    clusters_to_plot = range(n_gene_clusters - 1)

for cluster in clusters_to_plot:
    cluster_values = clustered_values[ordered_gene_clusters == cluster]
    n_cluster_genes = len(cluster_values)

    # Lay the axes out in inches: tissue strip on top, then one
    # CLUSTER_ROW_INCHES row per gene; labels sit outside and are kept by
    # bbox_inches="tight"
    strip_h = 0.3
    gap = 0.05
    heatmap_h = n_cluster_genes * CLUSTER_ROW_INCHES
    fig_h = heatmap_h + gap + strip_h
    fig = plt.figure(figsize=(CLUSTER_FIG_WIDTH, fig_h))

    ax_heatmap = fig.add_axes([0, 0, 1, heatmap_h / fig_h])
    ax_strip = fig.add_axes([0, (heatmap_h + gap) / fig_h, 1, strip_h / fig_h])

    ax_heatmap.imshow(
        cluster_values.to_numpy(),
        aspect="auto",
        interpolation="nearest",
        cmap=cmc.batlow,
        vmin=0,
        vmax=1
    )
    ax_heatmap.set_xticks(range(n_cols))
    ax_heatmap.set_xticklabels(cluster_values.columns, rotation=90)
    ax_heatmap.set_yticks(range(n_cluster_genes))
    ax_heatmap.set_yticklabels(
        gene_labels(cluster_values.index),
        fontsize=CLUSTER_GENE_FONT_SIZE
    )

    ax_strip.imshow(
        tissue_clusters[tissue_order][None, :],
        aspect="auto",
        interpolation="nearest",
        cmap=cluster_cmap(n_tissue_clusters),
        vmin=-0.5,
        vmax=n_tissue_clusters - 0.5
    )
    ax_strip.set_yticks([])
    ax_strip.set_xticks(range(n_cols))
    ax_strip.set_xticklabels(cluster_values.columns, rotation=90)
    ax_strip.xaxis.set_ticks_position("top")

    # Tissue names add about 2 in above and below the heatmap
    dpi = min(CLUSTER_DPI, int(MAX_PNG_PIXELS / (fig_h + 4)))
    file_stem = os.path.join(cluster_dir, f"PPI_ASIF_kmeans_cluster_{cluster + 1}")
    fig.savefig(f"{file_stem}.png", dpi=dpi, bbox_inches="tight")
    if fig_h + 4 <= MAX_PDF_INCHES:
        fig.savefig(f"{file_stem}.pdf", bbox_inches="tight")
    plt.close(fig)

    print(
        f"Cluster {cluster + 1}: {n_cluster_genes:,} genes, "
        f"{fig_h + 4:.0f} in tall, PNG at {dpi} dpi"
        + ("" if fig_h + 4 <= MAX_PDF_INCHES else ", no PDF")
    )