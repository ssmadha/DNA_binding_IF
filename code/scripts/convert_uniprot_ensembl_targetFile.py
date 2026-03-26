import pandas as pd
import biomart

TFLink_df = pd.read_csv("TFLink_Homo_sapiens_interactions_SS_simpleFormat_v1.0.tsv", sep="\t")
uniprot_ids = list(dict.fromkeys(TFLink_df.iloc[:,0].values.tolist() + TFLink_df.iloc[:,1].values.tolist()))

uniprot_ensg = {}
uniprot_symbol = {}

server = biomart.BiomartServer('http://useast.ensembl.org/biomart')
mart = server.datasets['hsapiens_gene_ensembl']
attributes = ['uniprot_gn_id', 'ensembl_gene_id', 'hgnc_symbol']
for i in range(0,len(uniprot_ids), 200):
    print(i)
    response = mart.search({'attributes': attributes,
                    'filters': {'uniprot_gn_id': uniprot_ids[i:min(i+200, len(uniprot_ids))]},
                   })
    data = response.raw.data.decode('ascii').strip().split('\n')
    for row in data:
        uniprot_id = row.split('\t')[0]
        ensg_id = row.split('\t')[1]
        symbol = row.split('\t')[2]
        uniprot_ensg[uniprot_id] = ensg_id
        uniprot_symbol[uniprot_id] = symbol

TFLink_df["TF"] = [uniprot_ensg[TF] if TF in uniprot_ensg else None for TF in TFLink_df.iloc[:,0]]
TFLink_df["target"] = [uniprot_ensg[target] if target in uniprot_ensg else None for target in TFLink_df.iloc[:,1]]
TFLink_df["target_symbol"] = [uniprot_symbol[target] if target in uniprot_ensg else None for target in TFLink_df.iloc[:,1]]

TFLink_df[["TF", "target", "target_symbol"]].to_csv("TFLink_Ensembl_ID.csv", index=False)