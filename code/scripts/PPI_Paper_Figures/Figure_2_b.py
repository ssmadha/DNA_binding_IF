import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


# ============================================================
# Load and prepare data
# ============================================================

df = pd.read_csv("ppi_binding_sites.tsv", sep="\t")

bc = df[df["Source"] == "bc"].copy()
w = df[df["Source"] == "w"].copy()
bc_w = df[df["Source"] == "bc_w"].copy()

# Add bc_w binding sites to BOTH bc and w
bc = pd.concat([bc, bc_w])
w = pd.concat([w, bc_w])

# Count binding sites per UniProt ID
bc_counts = bc.groupby("UniProt").size()
w_counts = w.groupby("UniProt").size()

# Total remains based on original data
total_counts = df.groupby("UniProt").size()


# ============================================================
# Venn diagram data
# ============================================================

bc_proteins = set(bc["UniProt"].dropna())
w_proteins = set(w["UniProt"].dropna())

bc_only = bc_proteins - w_proteins
w_only = w_proteins - bc_proteins
overlap = bc_proteins & w_proteins

n_bc_only = len(bc_only)
n_w_only = len(w_only)
n_overlap = len(overlap)
n_total = len(bc_proteins | w_proteins)


# ============================================================
# Histogram function
# ============================================================

def plot_binding_hist(ax, counts, title, fontsize=8):

    # Collapse everything above 20 into >20
    counts_plot = counts.clip(upper=21)

    # Bins for 1, 2, ..., 20, and >20
    bins = range(1, 23)

    ax.hist(
        counts_plot,
        bins=bins,
        align="left",
        rwidth=0.5       # Add space between bars
    )

    ax.set_title(title, fontsize=fontsize)

    ax.set_xlabel(
        "Number of Binding Sites",
        fontsize=fontsize
    )

    ax.set_ylabel(
        "Number of Genes",
        fontsize=fontsize
    )

    ax.tick_params(
        axis="both",
        labelsize=fontsize - 1
    )

    # Label the final bin as >20
    ticks = list(range(1, 21)) + [21]

    ax.set_xticks(ticks)
    ax.set_xticklabels(
        [str(i) for i in range(1, 21)] + [">20"]
    )

    ax.set_xlim(0.5, 21.5)

# ============================================================
# Main plot: Total
# ============================================================

fig, ax = plt.subplots(figsize=(12, 8))

plot_binding_hist(
    ax,
    total_counts,
    "Total Binding Sites per UniProt ID",
    fontsize=11
)


# ============================================================
# BC inset
# ============================================================

ax_bc = ax.inset_axes(
    [0.62, 0.62, 0.30, 0.32]
)

plot_binding_hist(
    ax_bc,
    bc_counts,
    "BC",
    fontsize=10
)


# ============================================================
# W inset
# ============================================================

ax_w = ax.inset_axes(
    [0.62, 0.20, 0.30, 0.32]
)

plot_binding_hist(
    ax_w,
    w_counts,
    "W",
    fontsize=10
)


# ============================================================
# Venn diagram
# ============================================================

ax_venn = ax.inset_axes(
    [0.32, 0.62, 0.30, 0.32]
)

ax_venn.set_xlim(0, 10)
ax_venn.set_ylim(0, 6)
ax_venn.axis("off")

bc_ellipse = Ellipse(
    (4.0, 3.2),
    width=5.0,
    height=4.0,
    alpha=0.4
)

w_ellipse = Ellipse(
    (6.0, 3.2),
    width=5.0,
    height=4.0,
    alpha=0.4
)

ax_venn.add_patch(bc_ellipse)
ax_venn.add_patch(w_ellipse)

ax_venn.text(
    2.8, 5.0,
    "BC",
    ha="center",
    va="center",
    fontsize=12,
    fontweight="bold"
)

ax_venn.text(
    7.2, 5.0,
    "W",
    ha="center",
    va="center",
    fontsize=12,
    fontweight="bold"
)

ax_venn.text(
    2.8, 3.2,
    str(n_bc_only),
    ha="center",
    va="center",
    fontsize=12
)

ax_venn.text(
    5.0, 3.2,
    str(n_overlap),
    ha="center",
    va="center",
    fontsize=12
)

ax_venn.text(
    7.2, 3.2,
    str(n_w_only),
    ha="center",
    va="center",
    fontsize=12
)

ax_venn.text(
    5.0, 0.5,
    f"Total unique proteins: {n_total}",
    ha="center",
    va="center",
    fontsize=10
)


plt.tight_layout()

plt.savefig("Figure_2_b.png")
plt.show()