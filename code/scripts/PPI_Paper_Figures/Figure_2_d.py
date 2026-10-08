import gzip
import os
import re

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import cmcrameri.cm as cmc
from matplotlib.colors import ListedColormap
from sklearn.cluster import KMeans

# ============================================================
# Settings
# ============================================================

NODE_DEGREE_CUTOFF = 50

# Second figure: genes with the highest mean ASIF across tissues,
# drawn from all genes regardless of node degree
TOP_ASIF_GENES = 50

# Third figure: all genes, k-means clustered (genes and tissues), with k
# chosen by the elbow method over these ranges
GENE_K_VALUES = range(3, 50)
TISSUE_K_VALUES = range(2, 21)
KMEANS_SEED = 0
# Leave the lowest-ASIF gene cluster (mostly ASIF ~ 0) out of the heatmap;
# clustering and the cluster tables still use all genes
DROP_LOWEST_GENE_CLUSTER = True
# Square heatmap tiles: side length in inches, and output resolution
TILE_INCHES = 0.1
KMEANS_DPI = 100


# ============================================================
# File paths
# ============================================================

ppi_file = "results_segment_ppi/expressed_coding_isoforms_with_relative_tpm_threshold_1_ASIF.tsv"
degree_file = "node_degree_df.tsv"
mapping_file = "mane_select_with_uniprot_id_mapping2.csv"
# Gene names (row labels) come from the same Ensembl GTF the pipeline uses
gtf_file = "reference_data/Homo_sapiens.GRCh38.109.gtf.gz"

output_dir = "ASIF_PPI_Figures"
os.makedirs(output_dir, exist_ok=True)


# ============================================================
# 1. Read files
# ============================================================

ppi_df = pd.read_csv(ppi_file, sep="\t")
degree_df = pd.read_csv(degree_file, sep="\t")
mapping_df = pd.read_csv(mapping_file)


# ============================================================
# 2. Map UniProt protein IDs -> Ensembl gene IDs
# ============================================================

mapping_filtered = mapping_df[
    ["gene_mane", "transcript_mane", "uniprot"]
].copy()

mapping_filtered = mapping_filtered.dropna(
    subset=["uniprot"]
)

degree_all_mapped = degree_df[
    ["protein", "node_degree"]
].merge(
    mapping_filtered,
    left_on="protein",
    right_on="uniprot",
    how="inner"
)

# Node degree per gene, for every mapped gene
gene_degree_all = (
    degree_all_mapped
    .groupby("gene_mane")["node_degree"]
    .max()
)


# ============================================================
# 3. Keep proteins above node-degree cutoff
# ============================================================

degree_mapped = degree_all_mapped[
    degree_all_mapped["node_degree"] >= NODE_DEGREE_CUTOFF
]

print(
    f"\nProteins with node degree >= {NODE_DEGREE_CUTOFF}: "
    f"{(degree_df['node_degree'] >= NODE_DEGREE_CUTOFF).sum()}"
)

print(
    f"Proteins successfully mapped to Ensembl: "
    f"{len(degree_mapped)}"
)


# ============================================================
# 4. Find ASIF columns
# ============================================================

asif_columns = [
    col for col in ppi_df.columns
    if col.endswith("_asif")
]

print(
    f"\nNumber of ASIF columns: "
    f"{len(asif_columns)}"
)


# ============================================================
# 5. Match Ensembl gene IDs to PPI table
# ============================================================

ppi_filtered = ppi_df[
    ppi_df["gene_id"].isin(
        degree_mapped["gene_mane"]
    )
].copy()

print(
    f"Rows in PPI table matching degree >= "
    f"{NODE_DEGREE_CUTOFF} proteins: "
    f"{len(ppi_filtered)}"
)


# ============================================================
# 6. Create heatmap data
# ============================================================

heatmap_df = ppi_filtered[
    ["gene_id"] + asif_columns
].copy()


# Take maximum ASIF across transcripts for each gene
heatmap_df = heatmap_df.groupby(
    "gene_id"
)[asif_columns].max()


# ============================================================
# 7. Clean tissue names
# ============================================================

heatmap_df.columns = [
    col.replace("_asif", "")
    for col in heatmap_df.columns
]


# ============================================================
# 8. Add node degree
# ============================================================

gene_degree = (
    degree_mapped
    .groupby("gene_mane")["node_degree"]
    .max()
)

heatmap_df["node_degree"] = (
    heatmap_df.index.map(gene_degree)
)

# Remove genes without a node-degree value
heatmap_df = heatmap_df.dropna(
    subset=["node_degree"]
)


# ============================================================
# 9. Sort by node degree, descending
# ============================================================

heatmap_df = heatmap_df.sort_values(
    "node_degree",
    ascending=False
)

node_degree = heatmap_df["node_degree"]

# ASIF values only
heatmap_values = heatmap_df.drop(
    columns="node_degree"
)


# ============================================================
# 9b. Gene names for row labels
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


def gene_labels(gene_ids):
    """Gene names for the given Ensembl IDs.

    Falls back to the Ensembl ID for genes without a name, and keeps the
    ID alongside any name shared by more than one gene in the figure.
    """
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


def plot_asif_heatmap(heatmap_values, node_degree, title, output_file):
    """ASIF heatmap (genes x tissues) with node-degree bars on the right.

    heatmap_values is indexed by Ensembl gene ID; node_degree is aligned
    to it and may be NaN for genes without a node degree (no bar drawn).
    """
    heatmap_values = heatmap_values.copy()
    heatmap_values.index = gene_labels(heatmap_values.index)

    # ========================================================
    # 10. Create figure
    # ========================================================

    n_genes = len(heatmap_values)

    fig = plt.figure(
        figsize=(18, max(8, n_genes * 0.25))
    )

    gs = fig.add_gridspec(
        nrows=1,
        ncols=2,
        width_ratios=[16, 3],
        wspace=0
    )

    ax_heatmap = fig.add_subplot(gs[0])
    ax_degree = fig.add_subplot(gs[1])

    # ========================================================
    # 11. ASIF heatmap
    # ========================================================

    sns.heatmap(
        heatmap_values,
        ax=ax_heatmap,
        cmap=cmc.batlow,
        vmin=0,
        vmax=1,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": "ASIF"}
    )

    ax_heatmap.set_xlabel("Tissue")
    ax_heatmap.set_ylabel("")

    ax_heatmap.set_title(title)

    ax_heatmap.tick_params(
        axis="x",
        rotation=90
    )

    # ========================================================
    # 12. Node-degree bars on the right
    # ========================================================

    y_positions = np.arange(n_genes) + 0.5

    ax_degree.barh(
        y_positions,
        node_degree.fillna(0).values,
        height=0.8
    )

    # Match heatmap vertical coordinates
    ax_degree.set_ylim(
        n_genes,
        0
    )

    # Remove vertical padding
    ax_degree.margins(y=0)

    # Remove y-axis
    ax_degree.set_yticks([])
    ax_degree.set_ylabel("")

    ax_degree.set_xlabel("Node degree")
    ax_degree.set_title("Node degree")

    # ========================================================
    # 13. Add node-degree values
    # ========================================================

    max_degree = node_degree.max()

    for y, value in zip(
        y_positions,
        node_degree.values
    ):

        ax_degree.text(
            0 if np.isnan(value) else value + max_degree * 0.015,
            y,
            "n/a" if np.isnan(value) else f"{int(value)}",
            va="center",
            ha="left",
            fontsize=8
        )

    # Give the labels some room
    ax_degree.set_xlim(
        0,
        max_degree * 1.15
    )

    # ========================================================
    # 14. Save
    # ========================================================

    plt.savefig(
        os.path.join(output_dir, output_file),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close(fig)


plot_asif_heatmap(
    heatmap_values,
    node_degree,
    title=(
        f"ASIF for Proteins with Node Degree ≥ "
        f"{NODE_DEGREE_CUTOFF}"
    ),
    output_file=(
        f"PPI_ASIF_node_degree_"
        f"{NODE_DEGREE_CUTOFF}_heatmap.png"
    )
)


# ============================================================
# 15. Genes with the highest mean ASIF across tissues
# ============================================================

# Same per-gene reduction as above (max over transcripts), for all genes
all_gene_asif = ppi_df.groupby(
    "gene_id"
)[asif_columns].max()

all_gene_asif.columns = [
    col.replace("_asif", "")
    for col in all_gene_asif.columns
]

mean_asif = all_gene_asif.mean(axis=1)
top_genes = mean_asif.sort_values(
    ascending=False
).index[:TOP_ASIF_GENES]

top_values = all_gene_asif.loc[top_genes]
top_degree = pd.Series(
    top_genes.map(gene_degree_all),
    index=top_genes
)

print(
    f"\nTop {TOP_ASIF_GENES} genes by mean ASIF: "
    f"mean ASIF {mean_asif[top_genes].min():.3f}-"
    f"{mean_asif[top_genes].max():.3f}, "
    f"{top_degree.isna().sum()} without a node degree"
)

plot_asif_heatmap(
    top_values,
    top_degree,
    title=(
        f"ASIF for the {TOP_ASIF_GENES} Genes with the "
        f"Highest Mean ASIF Across Tissues"
    ),
    output_file=f"PPI_ASIF_top_{TOP_ASIF_GENES}_mean_heatmap.png"
)


# ============================================================
# 16. All genes, k-means clustered (elbow method)
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

# Lay the axes out in inches so every heatmap tile is TILE_INCHES square
strip = 3 * TILE_INCHES
gap = TILE_INCHES
margin = 1.5
heatmap_w = n_cols * TILE_INCHES
heatmap_h = n_rows * TILE_INCHES
cbar_h = 0.12

fig_w = margin + strip + gap + heatmap_w + margin
fig_h = margin + heatmap_h + gap + strip + margin + cbar_h + 0.8

fig = plt.figure(figsize=(fig_w, fig_h))


def add_axes_inches(left, bottom, width, height):
    return fig.add_axes(
        [left / fig_w, bottom / fig_h, width / fig_w, height / fig_h]
    )


heatmap_left = margin + strip + gap
heatmap_bottom = margin
tissue_strip_bottom = heatmap_bottom + heatmap_h + gap

ax_heatmap = add_axes_inches(heatmap_left, heatmap_bottom, heatmap_w, heatmap_h)
ax_gene_strip = add_axes_inches(margin, heatmap_bottom, strip, heatmap_h)
ax_tissue_strip = add_axes_inches(heatmap_left, tissue_strip_bottom, heatmap_w, strip)
ax_cbar = add_axes_inches(
    heatmap_left, tissue_strip_bottom + strip + margin, heatmap_w, cbar_h
)

image = ax_heatmap.imshow(
    shown_values,
    aspect="auto",
    interpolation="nearest",
    cmap=cmc.batlow,
    vmin=0,
    vmax=1
)
ax_heatmap.set_xticks(range(n_cols))
ax_heatmap.set_xticklabels(clustered_values.columns, rotation=90, fontsize=7)
ax_heatmap.set_yticks([])
ax_heatmap.set_xlabel("Tissue")

colorbar = fig.colorbar(image, cax=ax_cbar, orientation="horizontal")
colorbar.ax.xaxis.set_ticks_position("top")
colorbar.ax.xaxis.set_label_position("top")
colorbar.set_label("ASIF")
colorbar.ax.tick_params(labelsize=7)

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
boundaries = np.flatnonzero(np.diff(shown_gene_clusters)) + 1
starts = np.concatenate(([0], boundaries))
ends = np.concatenate((boundaries, [n_rows]))
ax_gene_strip.set_yticks((starts + ends - 1) / 2)
ax_gene_strip.set_yticklabels(
    [str(cluster + 1) for cluster in shown_gene_clusters[starts]],
    fontsize=8
)
ax_gene_strip.set_ylabel(
    f"Genes (n = {n_rows:,}), k-means clusters"
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
ax_tissue_strip.set_xticklabels(clustered_values.columns, rotation=90, fontsize=7)
ax_tissue_strip.xaxis.set_ticks_position("top")
ax_tissue_strip.xaxis.set_label_position("top")

fig.suptitle(
    f"ASIF, k-means Clustered\n"
    f"({n_gene_clusters} gene / {n_tissue_clusters} tissue clusters)",
    x=(heatmap_left + heatmap_w / 2) / fig_w,
    y=1 - 0.1 / fig_h,
    va="top",
    fontsize=10
)

for extension in ["png", "pdf"]:
    fig.savefig(
        os.path.join(output_dir, f"PPI_ASIF_all_genes_kmeans_heatmap.{extension}"),
        dpi=KMEANS_DPI,
        bbox_inches="tight"
    )
plt.close(fig)

#plt.show()