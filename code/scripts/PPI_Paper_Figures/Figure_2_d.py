import gzip
import os
import re

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import cmcrameri.cm as cmc
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

# ============================================================
# Settings
# ============================================================

NODE_DEGREE_CUTOFF = 50

# Second figure: genes with the highest mean ASIF across tissues,
# drawn from all genes regardless of node degree
TOP_ASIF_GENES = 50

# Third figure: representative genes. "High" is ASIF >= HIGH_ASIF and
# "low" is ASIF < LOW_ASIF; each pattern contributes PER_PATTERN genes,
# chosen at random, and the rest of the RANDOM_GENES are drawn from all
# other genes
HIGH_ASIF = 0.8
LOW_ASIF = 0.2
FEW_TISSUES = (2, 5)
BROAD_MIN_TISSUES = 36
PER_PATTERN = 4
RANDOM_GENES = 8
REPRESENTATIVE_SEED = 42

FONT_SIZE = 14
plt.rcParams["font.size"] = FONT_SIZE
# Tissue names, and node-degree values and axis
TISSUE_FONT_SIZE = 17
DEGREE_FONT_SIZE = 16
# Shorter display names for some tissues
TISSUE_LABELS = {"parathyroid gland": "parathyroid"}


# ============================================================
# File paths
# ============================================================

ppi_file = "results_segment_ppi/expressed_coding_isoforms_with_relative_tpm_threshold_1_ASIF.tsv"
degree_file = "node_degree_df.tsv"
# Gene -> UniProt accessions, from the same Ensembl 109 xref file the
# pipeline uses (every protein of the gene, not only the MANE Select one)
mapping_file = "reference_data/Homo_sapiens.GRCh38.109.uniprot.tsv.gz"
# Gene names (row labels) come from the same Ensembl GTF the pipeline uses
gtf_file = "reference_data/Homo_sapiens.GRCh38.109.gtf.gz"

output_dir = "ASIF_PPI_Figures"
os.makedirs(output_dir, exist_ok=True)


# ============================================================
# 1. Read files
# ============================================================

ppi_df = pd.read_csv(ppi_file, sep="\t")
degree_df = pd.read_csv(degree_file, sep="\t")
mapping_df = pd.read_csv(mapping_file, sep="\t")


# ============================================================
# 2. Map UniProt protein IDs -> Ensembl gene IDs
# ============================================================

# One row per (gene, accession); an accession can belong to several
# genes (e.g. identical paralogs) and a gene to several accessions
mapping_filtered = (
    mapping_df[["gene_stable_id", "xref"]]
    .rename(columns={"gene_stable_id": "gene_id", "xref": "uniprot"})
    .drop_duplicates()
)

degree_all_mapped = degree_df[
    ["protein", "node_degree"]
].merge(
    mapping_filtered,
    left_on="protein",
    right_on="uniprot",
    how="inner"
)

# Node degree per gene (highest over its accessions), for every mapped gene
gene_degree_all = (
    degree_all_mapped
    .groupby("gene_id")["node_degree"]
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
    f"{degree_mapped['protein'].nunique()} "
    f"({degree_mapped['gene_id'].nunique()} genes)"
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
        degree_mapped["gene_id"]
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
    .groupby("gene_id")["node_degree"]
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


def plot_asif_heatmap(heatmap_values, node_degree, output_file):
    """ASIF heatmap (genes x tissues) with node-degree bars on the right.

    heatmap_values is indexed by Ensembl gene ID; node_degree is aligned
    to it and may be NaN for genes without a node degree (no bar drawn).
    The ASIF colorbar is saved separately as <output_file>_colorbar.png.
    """
    heatmap_values = heatmap_values.rename(columns=TISSUE_LABELS)
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
        wspace=0.02
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
        cbar=False
    )

    ax_heatmap.set_xlabel("")
    ax_heatmap.set_ylabel("")

    ax_heatmap.tick_params(
        axis="x",
        rotation=90,
        labelsize=TISSUE_FONT_SIZE
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
    ax_degree.tick_params(axis="x", labelsize=DEGREE_FONT_SIZE)

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
            fontsize=DEGREE_FONT_SIZE
        )

    # Give the labels some room
    ax_degree.set_xlim(
        0,
        max_degree * 1.35
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

    # ========================================================
    # 15. ASIF colorbar, saved on its own
    # ========================================================

    fig, ax_cbar = plt.subplots(figsize=(0.4, 4))
    colorbar = fig.colorbar(
        ScalarMappable(norm=Normalize(vmin=0, vmax=1), cmap=cmc.batlow),
        cax=ax_cbar
    )
    fig.savefig(
        os.path.join(
            output_dir,
            output_file.replace(".png", "_colorbar.png")
        ),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close(fig)


plot_asif_heatmap(
    heatmap_values,
    node_degree,
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
    output_file=f"PPI_ASIF_top_{TOP_ASIF_GENES}_mean_heatmap.png"
)


# ============================================================
# 16. Representative genes
# ============================================================

n_high = (all_gene_asif >= HIGH_ASIF).sum(axis=1)
n_low = (all_gene_asif < LOW_ASIF).sum(axis=1)
n_tissues = all_gene_asif.shape[1]

# High in the given number of tissues and low in all the others
patterns = {
    "one tissue": (n_high == 1) & (n_high + n_low == n_tissues),
    "few tissues": (
        n_high.between(*FEW_TISSUES) & (n_high + n_low == n_tissues)
    ),
    "across the board": n_high >= BROAD_MIN_TISSUES,
}

rng = np.random.default_rng(REPRESENTATIVE_SEED)
chosen = []
for pattern, is_pattern in patterns.items():
    candidates = all_gene_asif.index[is_pattern]
    picked = rng.choice(
        candidates,
        size=min(PER_PATTERN, len(candidates)),
        replace=False
    )
    print(f"Representative genes, {pattern}: {len(picked)} of {len(candidates)}")
    chosen += [(gene_id, pattern) for gene_id in picked]

already_chosen = {gene_id for gene_id, _ in chosen}
remaining = [g for g in all_gene_asif.index if g not in already_chosen]
picked = rng.choice(remaining, size=RANDOM_GENES, replace=False)
chosen += [(gene_id, "random") for gene_id in picked]

# Rows in order of increasing mean ASIF across tissues
representative = pd.DataFrame(chosen, columns=["gene_id", "pattern"])
representative["gene_name"] = representative["gene_id"].map(gene_names)
representative["n_tissues_high"] = representative["gene_id"].map(n_high)
representative["mean_asif"] = representative["gene_id"].map(mean_asif)
representative = representative.sort_values(
    "mean_asif", ascending=True, kind="stable"
)
representative["mean_asif"] = representative["mean_asif"].round(4)
representative.to_csv(
    os.path.join(output_dir, "PPI_ASIF_representative_genes.tsv"),
    sep="\t",
    index=False
)

representative_genes = pd.Index(representative["gene_id"])
plot_asif_heatmap(
    all_gene_asif.loc[representative_genes],
    pd.Series(
        representative_genes.map(gene_degree_all),
        index=representative_genes
    ),
    output_file="PPI_ASIF_representative_heatmap.png"
)

#plt.show()