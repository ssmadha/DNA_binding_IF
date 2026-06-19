import pandas as pd
import numpy as np

maradoner_df = pd.read_table("kulakovskiy_data/MARAdoner_activities.tsv", index_col=0)
maradoner_df = maradoner_df.set_index(maradoner_df.index.map(lambda x: x.split(".")[0]), append=True)
maradoner_df.index = maradoner_df.index.set_names(["MotifID", "GeneID"])

motif_info = pd.read_excel("kulakovskiy_data/SupplementaryTable11_MARA_Clusters.v2026.2.xlsx", sheet_name="TableS11")

asif_df = pd.read_table("ASIF_byGene_wName.csv", sep=",", index_col=[0,1])

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

maradoner_df_common = maradoner_df.loc[pd.IndexSlice[:, common_genes], common_columns.values()]

unique_motifs = motif_info.loc[motif_info["Cluster_Size"]==1, "Representative_Motif"].map(lambda x: x.split(".")[0]).tolist()

corrs = []
for val in common_columns.values():
    print(val)
    corr = asif_df_common.loc[pd.IndexSlice[:, unique_motifs], val].reset_index(drop=True).corr(maradoner_df_common[pd.IndexSlice[:, unique_motifs], val].reset_index(drop=True))
    corrs.append(asif_df_common[val].reset_index(drop=True).corr(maradoner_df_common[val].reset_index(drop=True)))

corrs = np.array(corrs)
corrs.mean()