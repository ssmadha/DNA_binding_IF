import os

import pandas as pd
import numpy as np

alpha = 63
beta = 0.3

#all_results.txt command:
#awk '{if($0~/ENSG/ && FNR!=1){printf "\n"$0} else if($0~/ENST/ || $0~/^[0-9]+$/){printf "\t"$1} else if($0~/^\[[0-9]/){printf "\t"$0}}' *.txt | awk '$3!=0' > all_results.txt
coverage_df = pd.read_csv("../../results_17k_redone/results/individual/all_results.txt", sep='\t', header=None)
coverage_df.columns = ["gene_id", "transcript_id", "n_domains", "domain_coverage"]
TPM_df = pd.read_csv("../../expressed_coding_isoforms_with_relative_tpm_threshold_1.tsv", sep='\t')

full_df = pd.merge(coverage_df, TPM_df)

def sigmoid(x):
    return 1/(1+np.power(np.e,-x))

asif_df = None
for i in range(full_df.shape[0]):
    curr_row = full_df.iloc[i, :]
    new_row = curr_row[:4].copy()
    coverages = [float(el) for el in curr_row["domain_coverage"][1:-1].split(", ")]
    impact_factor = 0
    for domain_coverage in coverages:
        # print(isoform_coverage_percentages[gene_id][transcript.protein][domain])
        # print(sigmoid(alpha*(1-isoform_coverage_percentages[gene_id][transcript.protein][domain]) \
        #                            + beta))
        impact_factor += sigmoid(alpha * (domain_coverage - beta))
    impact_factor /= curr_row["n_domains"]
    impact_factor = 1 - impact_factor
    for tpm in curr_row[4:]:
        tissue = curr_row.index[int(len(new_row)/2) + 2].split("_")[0]
        new_row = pd.concat([new_row, pd.Series({tissue+"_tpm": tpm, tissue+"_asif": impact_factor*tpm})])
    if asif_df is None:
        asif_df = pd.DataFrame([new_row])
    else:
        asif_df = pd.concat([asif_df, pd.DataFrame([new_row])])

asif_df.to_csv("../../PPI_ASIF_table2.tsv", sep='\t', index=False)