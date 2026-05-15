import pandas as pd

yue_df = pd.read_csv("../../merged_ppi.tsv", sep='\t')

contacts_a_df = yue_df[["protein_I_II", "protein_A", "contact_A", "source"]]
contacts_b_df = yue_df[["protein_I_II", "protein_B", "contact_B", "source"]]
contacts_a_df.columns = contacts_b_df.columns = ["ID", "UniProt", "Binding_Site", "Source"]

contacts_df = pd.concat([contacts_a_df, contacts_b_df], ignore_index=True).drop_duplicates()
contacts_df["ID"] = contacts_df["ID"] + "_" + contacts_df["UniProt"]

contacts_df.to_csv("../../ppi_binding_sites.tsv", sep='\t', index=False)