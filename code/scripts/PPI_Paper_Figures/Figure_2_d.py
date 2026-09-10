import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import cmcrameri.cm as cmc

# ============================================================
# Settings
# ============================================================

NODE_DEGREE_CUTOFF = 50


# ============================================================
# File paths
# ============================================================

ppi_file = "PPI_ASIF_table2.tsv"
degree_file = "node_degree_df.tsv"
mapping_file = "mane_select_with_uniprot_id_mapping2.csv"


# ============================================================
# 1. Read files
# ============================================================

ppi_df = pd.read_csv(ppi_file, sep="\t")
degree_df = pd.read_csv(degree_file, sep="\t")
mapping_df = pd.read_csv(mapping_file)


# ============================================================
# 2. Keep proteins above node-degree cutoff
# ============================================================

degree_filtered = degree_df.loc[
    degree_df["node_degree"] >= NODE_DEGREE_CUTOFF,
    ["protein", "node_degree"]
].copy()

print(
    f"\nProteins with node degree >= {NODE_DEGREE_CUTOFF}: "
    f"{len(degree_filtered)}"
)


# ============================================================
# 3. Map UniProt protein IDs -> Ensembl gene IDs
# ============================================================

mapping_filtered = mapping_df[
    ["gene_mane", "transcript_mane", "uniprot"]
].copy()

mapping_filtered = mapping_filtered.dropna(
    subset=["uniprot"]
)

degree_mapped = degree_filtered.merge(
    mapping_filtered,
    left_on="protein",
    right_on="uniprot",
    how="inner"
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
# 10. Create figure
# ============================================================

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


# ============================================================
# 11. ASIF heatmap
# ============================================================

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

ax_heatmap.set_title(
    f"ASIF for Proteins with Node Degree ≥ "
    f"{NODE_DEGREE_CUTOFF}"
)

ax_heatmap.tick_params(
    axis="x",
    rotation=90
)


# ============================================================
# 12. Node-degree bars on the right
# ============================================================

y_positions = np.arange(n_genes) + 0.5

ax_degree.barh(
    y_positions,
    node_degree.values,
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


# ============================================================
# 13. Add node-degree values
# ============================================================

max_degree = node_degree.max()

for y, value in zip(
    y_positions,
    node_degree.values
):

    ax_degree.text(
        value + max_degree * 0.015,
        y,
        f"{int(value)}",
        va="center",
        ha="left",
        fontsize=8
    )


# Give the labels some room
ax_degree.set_xlim(
    0,
    max_degree * 1.15
)


# ============================================================
# 14. Save
# ============================================================

output_file = (
    f"PPI_ASIF_node_degree_"
    f"{NODE_DEGREE_CUTOFF}_heatmap.png"
)

plt.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

#plt.show()