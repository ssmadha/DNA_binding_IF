import os

import matplotlib.pyplot as plt
import pandas as pd

print(os.getcwd())
target_df = pd.read_csv("../../visualization/asif/TFLink_Ensembl_ID.csv")
target_df.index = target_df["TF"]

asif_df = pd.read_csv("../../visualization/asif/ASIF_melted.csv")

joined_df = asif_df.join(target_df, on="Gene ID", how="left")


for tissue, df_sub in joined_df.groupby("Tissue"):

    df_sub = df_sub.groupby(["Gene ID", "ASIF"])[["Gene ID"]].apply(lambda x: len(x)).reset_index()
    print(df_sub.head())
    plt.figure()
    plt.plot(df_sub.iloc[0], df_sub["ASIF"], ".")
    plt.title(f"Tissue: {tissue}")
    plt.xlabel("Number of targets")
    plt.ylabel("ASIF")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()
