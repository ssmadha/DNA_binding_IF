import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle

maradoner_df = pd.read_table("kulakovskiy_data/MARAdoner_activities.tsv", index_col=0)
maradoner_df = maradoner_df.set_index(maradoner_df.index.map(lambda x: x.split(".")[0]), append=True)
maradoner_df.index = maradoner_df.index.set_names(["MotifID", "GeneID"])

motif_info = pd.read_excel("kulakovskiy_data/SupplementaryTable11_MARA_Clusters.v2026.2.xlsx", sheet_name="TableS11")

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

asif_transcript = pickle.load(open("dbd_impact_factors1.pickle", "rb"))

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

asif_df_common = asif_df.loc[pd.IndexSlice[:, common_genes], common_columns.keys()]
asif_df_common = asif_df_common.rename(columns=common_columns)

asif_transcript_df_common = asif_transcript_df.loc[:, common_columns.keys()]

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
for gene in common_genes:
    mask_asif = asif_df_common.index.get_level_values("Gene Name").isin([gene])
    mask_maradoner = maradoner_df_common.index.get_level_values("GeneID").isin([gene])

    print(gene)

    x = asif_df_common.loc[mask_asif, :].to_numpy().ravel()
    y = maradoner_df_common.loc[mask_maradoner, :].mean(axis=0).to_numpy()

    # Remove NaN pairs
    valid = ~(pd.isna(x) | pd.isna(y))
    x = x[valid]
    y = y[valid]

    corr = pd.Series(x).corr(pd.Series(y))
    corrs[gene] = corr

    print(f"{gene}: {corr:.3f}")

    plt.figure(figsize=(5, 5))
    plt.scatter(x, y)
    gene_id = gene_name_to_id.get(gene)

    if gene_id is not None:
        mask_transcripts = (
                asif_transcript_df_common.index.get_level_values("Gene ID") == gene_id
        )

        transcript_df = asif_transcript_df_common.loc[mask_transcripts]

        for transcript_id, row in transcript_df.iterrows():
            xt = row.to_numpy()

            valid = ~(pd.isna(xt) | pd.isna(y))

            plt.scatter(
                xt[valid],
                y[valid],
                s=20,
                alpha=0.6,
                marker="x",
                label="Transcript" if transcript_id == transcript_df.index[0] else None,
            )

    plt.xlim([0, 1])
    plt.xlabel("ASIF")
    plt.ylabel("Maradoner (mean)")
    plt.title(f"{gene}\nPearson r = {corr:.3f}; x.var() = {x.var()}")
    plt.tight_layout()

    save_folder = ""
    if x.var() < 0.001:
        save_folder = "Category6/"
        corrs6[gene] = corr
    #Category 1: As expected negatively correlated
    elif corr < -0.1:
        save_folder = "Category1/"
        corrs1[gene] = corr
    #Category 2: Anti-expected positive correlated
    elif corr > 0.1:
        save_folder = "Category2/"
        corrs2[gene] = corr
    #Category 3: Not much change at 0
    elif -3 < y.mean() < 3:
        save_folder = "Category3/"
        corrs3[gene] = corr
    #Category 4: Not much change positive
    elif y.mean() > 3:
        save_folder = "Category4/"
        corrs4[gene] = corr
    #Category 5: Not much change negative
    elif y.mean() < -3:
        save_folder = "Category5/"
        corrs5[gene] = corr
    #Category 6: random
    else:
        save_folder = "Category6/"
        corrs6[gene] = corr

    # Save plot
    plt.savefig(f"kulakovskiy_data/corr_graphs/{save_folder}{gene}_scatter.png", dpi=300, bbox_inches="tight")

    plt.close()


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
