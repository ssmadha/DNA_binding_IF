import pandas as pd
import matplotlib.pyplot as plt

# Load the TSV
df = pd.read_csv("ppi_binding_sites.tsv", sep="\t")

# Get the three sources
sources = df["Source"].dropna().unique()

# Create figure with 4 histograms
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()

# Plot each source
for ax, source in zip(axes[:3], sources):

    # Count binding sites per UniProt ID for this source
    binding_sites = (
        df[df["Source"] == source]
        .groupby("UniProt")
        .size()
    )

    # Histogram
    ax.hist(
        binding_sites,
        bins=range(1, binding_sites.max() + 2),
        align="left"
    )

    ax.set_title(source)
    ax.set_xlabel("Number of Binding Sites")
    ax.set_ylabel("Number of Genes")

# Total across all sources
total_binding_sites = (
    df.groupby("UniProt")
      .size()
)

axes[3].hist(
    total_binding_sites,
    bins=range(1, total_binding_sites.max() + 2),
    align="left"
)

axes[3].set_title("Total")
axes[3].set_xlabel("Number of Binding Sites")
axes[3].set_ylabel("Number of Genes")

plt.tight_layout()
plt.show()