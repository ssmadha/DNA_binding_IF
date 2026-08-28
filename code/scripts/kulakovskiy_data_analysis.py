import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.stats import pearsonr, gaussian_kde, binomtest, kruskal, mannwhitneyu, chi2_contingency
from statsmodels.stats.multitest import multipletests

#don't automatically show the plots
plt.ioff()

maradoner_df = pd.read_table("kulakovskiy_data/MARAdoner_activities.tsv", index_col=0)
maradoner_df = maradoner_df.set_index(maradoner_df.index.map(lambda x: x.split(".")[0]), append=True)
maradoner_df.index = maradoner_df.index.set_names(["MotifID", "GeneID"])

motif_info = pd.read_excel("kulakovskiy_data/SupplementaryTable11_MARA_Clusters.v2026.2.xlsx", sheet_name="TableS11")

motif_families = dict()
for i in motif_info.index:
    this_rep_motif = motif_info.loc[i, "Representative_Motif"]
    family_members = motif_info.loc[i, "Clustered_Motifs_Gene_Symbols"].split("; ")
    motif_families[this_rep_motif] = family_members


asif_df = pd.read_table("visualization/asif/ASIF.csv", sep=",", index_col=[0, 1])

asif_DBD = pickle.load(open("dbd_impact_factors1.pickle", "rb"))
asif_DBD_df = pd.DataFrame.from_dict(asif_DBD, orient="index")
gene_map = (
    asif_df.index.to_frame(index=False)
    .drop_duplicates(subset="Gene ID")
    .set_index("Gene ID")["Gene Name"]
)

# Add Gene ID column to the new dataframe
asif_DBD_df["Gene Name"] = asif_DBD_df.index.map(gene_map)

# Make a MultiIndex
asif_DBD_df = (
    asif_DBD_df
    .reset_index(names="Gene ID")
    .set_index(["Gene ID", "Gene Name"])
)

asif_transcript = pickle.load(open("rel_impact_factors1.pickle", "rb"))

asif_rows = []
exp_rows = []

for gene_id, transcripts in asif_transcript.items():
    for transcript_id, tissues in transcripts.items():
        asif_row = {
            "Gene ID": gene_id,
            "Transcript ID": transcript_id
        }
        exp_row = {
            "Gene ID": gene_id,
            "Transcript ID": transcript_id
        }

        for tissue, (asif_value, exp_value) in tissues.items():
            asif_row[tissue] = asif_value
            exp_row[tissue] = exp_value

        asif_rows.append(asif_row)
        exp_rows.append(exp_row)

asif_transcript_df = (
    pd.DataFrame(asif_rows)
      .set_index(["Gene ID", "Transcript ID"])
      .sort_index()
)

exp_transcript_df = (
    pd.DataFrame(exp_rows)
      .set_index(["Gene ID", "Transcript ID"])
      .sort_index()
)

# Add Gene ID column to the new dataframe
asif_transcript_df["Gene Name"] = asif_transcript_df.index.get_level_values("Gene ID").map(gene_map)

# Make a MultiIndex
asif_transcript_df = (
    asif_transcript_df
    .reset_index(names=["Gene ID", "Transcript ID"])
    .set_index(["Gene ID", "Transcript ID", "Gene Name"])
)

common_columns = {"breast": "breast_adult",
                  "epididymis": "epididymis_adult",
                  "appendix": "appendix_adult",
                  "spleen": "spleen_adult",
                  "skin": "skin_adult",
                  "tonsil": "tonsil_adult",
                  "thymus": "thymus_adult",
                  "salivary gland": "salivary_gland_adult",
                  "prostate": "prostate_adult",
                  "kidney": "kidney_adult",
                  "adipose tissue": "adipose_tissue_adult",
                  "urinary bladder": "bladder_adult",
                  "esophagus": "esophagus_adult",
                  "thyroid gland": "thyroid_adult",
                  "adrenal gland": "adrenal_gland_adult",
                  "smooth muscle": "smooth_muscle_adult",
                  "cervix": "cervix_adult",
                  "gallbladder": "gall_bladder_adult",
                  "colon": "colon_adult",
                  "seminal vesicle": "seminal_vesicle_adult",
                  "placenta": "placenta_adult",
                  "ovary": "ovary_adult",
                  "lung": "lung_adult",
                  "liver": "liver_adult",
                  "skeletal muscle": "skeletal_muscle_adult",
                  "small intestine": "small_intestine_adult",
                  "pancreas": "pancreas_adult",
                  "testis": "testis_adult",
                  "bone marrow": "bone_marrow_adult",
                  "tongue": "tongue_adult"
}

common_genes = maradoner_df.index.get_level_values(1).intersection(asif_df.index.get_level_values(1)).tolist()

asif_df_common = asif_df.loc[:, common_columns.keys()]
asif_df_common = asif_df_common.rename(columns=common_columns)

asif_transcript_df_common = asif_transcript_df.loc[:, common_columns.keys()]
asif_transcript_df_common = (
    asif_transcript_df
    .loc[:, common_columns.keys()]
    .rename(columns=common_columns)
)

maradoner_df_common = maradoner_df.loc[pd.IndexSlice[:, common_genes], common_columns.values()]

gene_name_to_id = (
    asif_df.index.to_frame(index=False)
    .drop_duplicates("Gene Name")
    .set_index("Gene Name")["Gene ID"]
)

corrs = dict()
corrs1 = dict()
corrs2 = dict()
corrs3 = dict()
corrs4 = dict()
corrs5 = dict()
corrs6 = dict()
for rep_motif, family_members in motif_families.items():
    print(rep_motif)
    # Skip motifs that aren't present in the Maradoner dataframe
    mask_maradoner = (
        maradoner_df_common.index.get_level_values("MotifID") == rep_motif
    )

    if not mask_maradoner.any():
        continue

    y = maradoner_df_common.loc[mask_maradoner].mean(axis=0).to_numpy()

    max_x = float("-inf")

    plt.figure(figsize=(6, 6))

    for family_gene in family_members:
        print(family_gene)
        mask_asif = (
            asif_df_common.index.get_level_values("Gene Name") == family_gene
        )

        if not mask_asif.any():
            continue

        x = asif_df_common.loc[mask_asif].mean(axis=0).to_numpy()

        valid = ~(np.isnan(x) | np.isnan(y))

        x_valid = x[valid]
        y_valid = y[valid]

        if len(x_valid) < 2:
            continue

        corr, pval = pearsonr(x_valid, y_valid)
        corrs[(rep_motif, family_gene)] = {
            "r": corr,
            "p": pval
        }

        max_x = max(max_x, np.max(x_valid))

        plt.scatter(
            x_valid,
            y_valid,
            label=f"{family_gene} (r={corr:.2f}, p={pval:.2g})"
        )

        gene_id = gene_name_to_id.get(family_gene)

        # if gene_id is not None:
        #
            # mask_transcripts = (
            #         asif_transcript_df_common.index.get_level_values("Gene ID")
            #         == gene_id
            # )
            #
            # transcript_df = asif_transcript_df_common.loc[mask_transcripts]
            #
            # for (_, transcript_id, _), row in transcript_df.iterrows():
            #     xt = row.to_numpy()
            #
            #     valid = ~(np.isnan(xt) | np.isnan(y))
            #
            #     plt.scatter(
            #         xt[valid],
            #         y[valid],
            #         marker="x",
            #         s=20,
            #         alpha=0.5,
            #         label=transcript_id
            #     )

    plt.xlabel("ASIF")
    plt.ylabel("Maradoner")
    plt.title(rep_motif)
    plt.legend(
        fontsize=7,
        bbox_to_anchor=(1.02, 1),
        loc="upper left"
    )

    label_x = max_x + 0.03  # adjust spacing as needed

    for tissue, yi in zip(common_columns.keys(), y):
        if np.isnan(yi):
            continue

        plt.text(
            label_x,
            yi,
            tissue,
            fontsize=7,
            va="center",
            ha="left"
        )
    xmin, xmax = plt.xlim()
    plt.xlim(xmin, label_x + 0.05)

    plt.tight_layout()

    save_folder = ""
    if x.var() < 0.01:
        save_folder = "Category6/"
        corrs6[rep_motif] = corr
    #Category 1: As expected negatively correlated
    elif corr < -0.1:
        save_folder = "Category1/"
        corrs1[rep_motif] = corr
    #Category 2: Anti-expected positive correlated
    elif corr > 0.1:
        save_folder = "Category2/"
        corrs2[rep_motif] = corr
    #Category 3: Not much change at 0
    elif -2 < y.mean() < 2:
        save_folder = "Category3/"
        corrs3[rep_motif] = corr
    #Category 4: Not much change positive
    elif y.mean() > 2:
        save_folder = "Category4/"
        corrs4[rep_motif] = corr
    #Category 5: Not much change negative
    elif y.mean() < -2:
        save_folder = "Category5/"
        corrs5[rep_motif] = corr
    #Category 6: random
    else:
        save_folder = "Category6/"
        corrs6[rep_motif] = corr

    # Save plot
    plt.savefig(f"kulakovskiy_data/corr_graphs_bymotif/{save_folder}{rep_motif}_scatter.png", dpi=300, bbox_inches="tight")

    plt.savefig(f"kulakovskiy_data/corr_graphs_bymotif/{rep_motif}_scatter.png", dpi=300, bbox_inches="tight")

    plt.close()

plt.figure(figsize=(6, 6))
for rep_motif, family_members in motif_families.items():
    #print(rep_motif)
    # Skip motifs that aren't present in the Maradoner dataframe
    mask_maradoner = (
        maradoner_df_common.index.get_level_values("MotifID") == rep_motif
    )

    if not mask_maradoner.any():
        continue

    y = maradoner_df_common.loc[mask_maradoner].mean(axis=0).to_numpy()

    # if -2 < y.mean() < 2:
    #     continue

    for family_gene in family_members:
        #print(family_gene)
        mask_asif = (
            asif_df_common.index.get_level_values("Gene Name") == family_gene
        )

        if not mask_asif.any():
            continue

        x = asif_df_common.loc[mask_asif].mean(axis=0).to_numpy()

        valid = ~(np.isnan(x) | np.isnan(y))

        x_valid = x[valid]
        y_valid = y[valid]

        no_zero = (-2 > y_valid) | (y_valid > 2)

        x_no_zero = x_valid[no_zero]
        y_no_zero = y_valid[no_zero]

    if len(x_no_zero) < 2:
        continue

    plt.scatter(
        x_valid,
        y_valid,
        label=f"{family_gene}"
    )

plt.xlabel("ASIF")
plt.ylabel("MARADONER")
plt.tight_layout()
plt.show()

plt.close()

all_x = []
all_y = []

for rep_motif, family_members in motif_families.items():

    mask_maradoner = (
        maradoner_df_common.index
        .get_level_values("MotifID") == rep_motif
    )

    if not mask_maradoner.any():
        continue

    y = maradoner_df_common.loc[mask_maradoner].mean(axis=0).to_numpy()

    if -2 < y.mean() < 2:
        continue

    for family_gene in family_members:

        mask_asif = (
            asif_df_common.index
            .get_level_values("Gene Name") == family_gene
        )

        if not mask_asif.any():
            continue

        x = asif_df_common.loc[mask_asif].mean(axis=0).to_numpy()

        valid = ~(np.isnan(x) | np.isnan(y))

        x_valid = x[valid]
        y_valid = y[valid]

        no_zero = (-2 > y_valid) | (y_valid > 2)

        x_no_zero = x_valid[no_zero]
        y_no_zero = y_valid[no_zero]

        all_x.extend(x_no_zero)
        all_y.extend(y_no_zero)

# Convert to arrays
all_x = np.asarray(all_x)
all_y = np.asarray(all_y)

# KDE
xy = np.vstack([all_x, all_y])
kde = gaussian_kde(xy)

# Create grid
xmin, xmax = all_x.min(), all_x.max()
ymin, ymax = all_y.min(), all_y.max()

#ymin, ymax = 10, -10

xx, yy = np.meshgrid(
    np.linspace(xmin, xmax, 200),
    np.linspace(ymin, ymax, 200)
)

grid_coords = np.vstack([xx.ravel(), yy.ravel()])
density = kde(grid_coords).reshape(xx.shape)

levels = np.linspace(
    density.min(),
    density.max(),
    40
)

plt.figure(figsize=(7, 7))

contours = plt.contour(
    xx,
    yy,
    density,
    levels=levels
)

plt.clabel(contours, inline=True, fontsize=7)

plt.xlabel("ASIF")
plt.ylabel("MARADONER")
plt.title("ASIF vs MARADONER point density")

plt.tight_layout()
plt.show()

rows = []

for rep_motif, family_members in motif_families.items():

    # Maradoner rows for this representative motif
    mask_maradoner = (
        maradoner_df_common.index.get_level_values("MotifID") == rep_motif
    )

    if not mask_maradoner.any():
        continue

    maradoner_subset = maradoner_df_common.loc[mask_maradoner]

    # Mean Maradoner value across rows for each tissue
    maradoner_mean = maradoner_subset.mean(axis=0)

    for family_gene in family_members:

        # ASIF gene-level data
        mask_asif = (
            asif_df_common.index.get_level_values("Gene Name")
            == family_gene
        )

        if not mask_asif.any():
            continue

        asif_subset = asif_df_common.loc[mask_asif]

        # If there are multiple ASIF rows for the gene, take the mean
        asif_values = asif_subset.mean(axis=0)

        # Gene ID
        gene_id = gene_name_to_id.get(family_gene)

        # Transcript-level ASIF
        if gene_id is not None:

            mask_transcripts = (
                asif_transcript_df_common.index
                .get_level_values("Gene ID") == gene_id
            )

            transcript_subset = asif_transcript_df_common.loc[
                mask_transcripts
            ]

        else:
            transcript_subset = pd.DataFrame()

        # One row for every tissue
        for tissue in common_columns.values():

            asif_value = asif_values.get(tissue, np.nan)
            maradoner_value = maradoner_mean.get(tissue, np.nan)

            # Find transcript with maximum ASIF for this tissue
            if not transcript_subset.empty and tissue in transcript_subset.columns:

                transcript_values = transcript_subset[tissue]

                if transcript_values.notna().any():

                    max_transcript_idx = transcript_values.idxmax()
                    max_asif = transcript_values.loc[max_transcript_idx]

                    # MultiIndex -> extract Transcript ID
                    if isinstance(max_transcript_idx, tuple):
                        transcript_id = max_transcript_idx[1]
                    else:
                        transcript_id = max_transcript_idx

                    # Use transcript maximum rather than gene-level value
                    asif_max = max_asif

                else:
                    asif_max = np.nan
                    transcript_id = np.nan

            else:
                asif_max = asif_value
                transcript_id = np.nan

            rows.append({
                "Maradoner": maradoner_value,
                "ASIF_MAX": asif_max,
                "Gene name": family_gene,
                "Gene ID": gene_id,
                "Tissue": tissue,
                "Specific Isoform that has the Max ASIF": transcript_id
            })

combined_df = pd.DataFrame(rows)

with open("DNA_binding_lost.pickle", "rb") as handle:
    DNA_binding_lost = pickle.load(handle)

combined_df["DNA_binding_lost"] = combined_df.apply(
    lambda row: DNA_binding_lost
        .get(row["Gene ID"], {})
        .get(row["Specific Isoform that has the Max ASIF"], np.nan),
    axis=1
)

combined_df.to_csv("Maradoner_ASIF_table.tsv", sep="\t", index=False)

subset_near_zero = combined_df[
    combined_df["Maradoner"].between(-2, 2, inclusive="neither")
].copy()

subset_near_zero["ASIF_category"] = np.where(
    subset_near_zero["ASIF_MAX"] >= 0.5,
    "High ASIF",
    "Low ASIF"
)

counts = subset_near_zero["ASIF_category"].value_counts()

high = counts.get("High ASIF", 0)
low = counts.get("Low ASIF", 0)

result = binomtest(
    high,
    n=high + low,
    p=0.5,
    alternative="two-sided"
)

print(f"High ASIF: {high}")
print(f"Low ASIF: {low}")
print(f"p-value: {result.pvalue:.4g}")

subset_far_zero = combined_df[
    ~combined_df["Maradoner"].between(-2, 2, )
].copy()

subset_far_zero["ASIF_category"] = np.where(
    subset_far_zero["ASIF_MAX"] >= 0.5,
    "High ASIF",
    "Low ASIF"
)

print(subset_far_zero["ASIF_category"].value_counts())

subset = combined_df.dropna(subset=["Maradoner", "ASIF_MAX"]).copy()

subset["Maradoner_bin"] = pd.cut(
    subset["Maradoner"],
    bins=[-np.inf, -5, -2, 0, 2, 5, np.inf],
    labels=["< -5", "-5 to -2", "-2 to 0", "0 to 2", "2 to 5", "> 5"]
)

summary = (
    subset.groupby("Maradoner_bin", observed=True)["ASIF_MAX"]
    .agg(
        count="count",
        mean="mean",
        median="median",
        std="std",
        below_0_2=lambda x: (x < 0.2).sum(),
        above_0_8=lambda x: (x >= 0.8).sum()
    )
)
summary["pct_below_0_2"] = (
    summary["below_0_2"] / summary["count"] * 100
)

summary["pct_above_0_8"] = (
    summary["above_0_8"] / summary["count"] * 100
)

print(summary)

bin_order = [
    "< -5",
    "-5 to -2",
    "-2 to 0",
    "0 to 2",
    "2 to 5",
    "> 5"
]

# Get ASIF values for each bin
box_data = [
    subset.loc[
        subset["Maradoner_bin"] == bin_name,
        "ASIF_MAX"
    ].values
    for bin_name in bin_order
]

plt.figure(figsize=(9, 6))

plt.boxplot(
    box_data,
    tick_labels=bin_order,
    showfliers=False
)

# Add jittered individual observations
for i, values in enumerate(box_data, start=1):
    jitter = np.random.normal(
        i,
        0.05,
        size=len(values)
    )

    plt.scatter(
        jitter,
        values,
        alpha=0.25,
        s=12
    )

plt.ylim(0, 1)

plt.xlabel("Maradoner score")
plt.ylabel("ASIF")
plt.title("ASIF distribution across Maradoner score ranges")

plt.tight_layout()
plt.show()

# Get ASIF values for each Maradoner bin
violin_data = [
    subset.loc[
        subset["Maradoner_bin"] == bin_name,
        "ASIF_MAX"
    ].values
    for bin_name in bin_order
]

plt.figure(figsize=(9, 6))

parts = plt.violinplot(
    violin_data,
    positions=np.arange(1, len(bin_order) + 1),
    showmeans=False,
    showmedians=True,
    showextrema=True
)

# ASIF ranges from 0 to 1
plt.ylim(0, 1)

# Reference line at ASIF = 0.5
plt.axhline(
    0.5,
    linestyle="--",
    linewidth=1
)

plt.xticks(
    np.arange(1, len(bin_order) + 1),
    bin_order
)

plt.xlabel("Maradoner score")
plt.ylabel("ASIF")
plt.title("ASIF distribution across Maradoner score ranges")

plt.tight_layout()
plt.show()

# Data for each group
groups = {
    b: subset.loc[
        subset["Maradoner_bin"] == b,
        "ASIF_MAX"
    ].values
    for b in bin_order
}

# Remove empty groups
groups = {
    k: v for k, v in groups.items()
    if len(v) > 0
}

# --------------------------------------------------
# Overall Kruskal-Wallis test
# --------------------------------------------------

kw_stat, kw_p = kruskal(*groups.values())

print("Kruskal-Wallis test")
print(f"H = {kw_stat:.3f}")
print(f"p = {kw_p:.4g}")


# --------------------------------------------------
# Pairwise Mann-Whitney U tests
# --------------------------------------------------

results = []

group_names = list(groups.keys())

for i in range(len(group_names)):
    for j in range(i + 1, len(group_names)):

        group1 = group_names[i]
        group2 = group_names[j]

        x = groups[group1]
        y = groups[group2]

        stat, p = mannwhitneyu(
            x,
            y,
            alternative="two-sided"
        )

        results.append({
            "Group 1": group1,
            "Group 2": group2,
            "n1": len(x),
            "n2": len(y),
            "U": stat,
            "p": p
        })

results_df = pd.DataFrame(results)

# --------------------------------------------------
# FDR correction
# --------------------------------------------------

results_df["p_adj"] = multipletests(
    results_df["p"],
    method="fdr_bh"
)[1]

results_df["Significant"] = results_df["p_adj"] < 0.05

print(results_df)

# --------------------------------------------------
# Chi^2 test
# --------------------------------------------------

# Categorize ASIF
subset["ASIF_category"] = pd.cut(
    subset["ASIF_MAX"],
    bins=[-np.inf, 0.2, 0.8, np.inf],
    labels=["Low (<0.2)", "Intermediate (0.2-0.8)", "High (>=0.8)"]
)

# Contingency table
contingency = pd.crosstab(
    subset["Maradoner_bin"],
    subset["ASIF_category"]
)

print(contingency)

chi2, p, dof, expected = chi2_contingency(contingency)

print(f"Chi-square = {chi2:.3f}")
print(f"Degrees of freedom = {dof}")
print(f"p-value = {p:.4g}")

print(np.array(list(corrs.values())).mean())

unique_motifs = motif_info.loc[motif_info["Cluster_Size"]<2, "Representative_Motif"].map(lambda x: x.split(".")[0]).tolist()

unique_corrs = {str(gene): corrs[gene] for gene in unique_motifs if gene in common_genes}

print(np.array(list(unique_corrs.values())).mean())


asif_DBD_df_common = asif_df.loc[pd.IndexSlice[:, common_genes], common_columns.keys()]
asif_DBD_df_common = asif_DBD_df_common.rename(columns=common_columns)

corrs_DBD = dict()
for gene in common_genes:
    mask_asif = asif_DBD_df_common.index.get_level_values("Gene Name").isin([gene])
    mask_maradoner = maradoner_df_common.index.get_level_values("GeneID").isin([gene])
    print(gene)
    x = asif_DBD_df_common.loc[mask_asif, :].to_numpy().ravel()
    y = maradoner_df_common.loc[mask_maradoner, :].mean(axis=0).to_numpy()

    # Remove NaN pairs
    valid = ~(pd.isna(x) | pd.isna(y))
    x = x[valid]
    y = y[valid]

    corr = pd.Series(x).corr(pd.Series(y))
    corrs_DBD[gene] = corr

    print(f"{gene}: {corr:.3f}")

    plt.figure(figsize=(5, 5))
    plt.scatter(x, y)
    plt.xlabel("ASIF")
    plt.ylabel("Maradoner (mean)")
    plt.title(f"{gene}\nPearson r = {corr:.3f}")
    plt.tight_layout()

    # Save plot
    plt.savefig(f"kulakovskiy_data/corr_graphs_DBD/{gene}_scatter.png", dpi=300, bbox_inches="tight")

    plt.close()


print(np.array(list(corrs_DBD.values())).mean())

unique_motifs = motif_info.loc[motif_info["Cluster_Size"]<2, "Representative_Motif"].map(lambda x: x.split(".")[0]).tolist()

unique_corrs_DBD = {str(gene): corrs_DBD[gene] for gene in unique_motifs if gene in common_genes}

print(np.array(list(unique_corrs_DBD.values())).mean())
