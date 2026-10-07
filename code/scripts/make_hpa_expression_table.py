#!/usr/bin/env python
"""
Build a transcript expression table for compute_asif.py from the Human
Protein Atlas transcript tissue RNA file (transcript_rna_tissue.tsv).

The HPA file has ensgid, enstid, then per-sample "TPM.<tissue>.<n>" and
"est_counts.<tissue>.<n>" columns. Samples are averaged per tissue (as in
getTFIsoforms.ipynb), and with the default scaling each transcript's mean
TPM is divided by the highest mean TPM among its gene's transcripts in that
tissue (all HPA transcripts of the gene, so the top isoform is 1). Tissues
where the gene has no expression get 0.

Output: TSV with gene_id, transcript_id, then "<tissue>_TPM" per tissue in
the HPA file's column order.
"""
import argparse

import pandas as pd


def tissue_mean_tpm(hpa_file):
    """
    Read the HPA file and average its TPM samples per tissue.

    Parameters
    ----------
    hpa_file: str
        Path to transcript_rna_tissue.tsv.

    Returns
    -------
    pandas.DataFrame
        gene_id, transcript_id, then one mean-TPM column per tissue.
    """
    header = pd.read_csv(hpa_file, sep="\t", nrows=0).columns
    tpm_columns = [column for column in header if column.startswith("TPM.")]
    hpa = pd.read_csv(hpa_file, sep="\t", usecols=["ensgid", "enstid"] + tpm_columns)
    tissues = list(dict.fromkeys(column.split(".")[1] for column in tpm_columns))
    means = pd.DataFrame({tissue: hpa[[column for column in tpm_columns
                                       if column.split(".")[1] == tissue]].mean(axis=1)
                          for tissue in tissues})
    means.insert(0, "transcript_id", hpa["enstid"])
    means.insert(0, "gene_id", hpa["ensgid"])
    return means


def relative_to_max(means):
    """
    Divide each transcript's per-tissue TPM by its gene's highest
    transcript TPM in that tissue; 0 where the gene's maximum is 0.

    Parameters
    ----------
    means: pandas.DataFrame
        Output of tissue_mean_tpm.

    Returns
    -------
    pandas.DataFrame
    """
    tissues = means.columns[2:]
    gene_max = means.groupby("gene_id")[tissues].transform("max")
    relative = means.copy()
    relative[tissues] = (means[tissues] / gene_max).fillna(0)
    return relative


def get_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", default="../transcript_rna_tissue.tsv",
                        help="HPA transcript_rna_tissue.tsv. Default: %(default)s")
    parser.add_argument("-o", "--output", required=True,
                        help="Path to write the expression TSV to")
    parser.add_argument("-s", "--scaling", choices=["relative_to_max", "absolute"],
                        default="relative_to_max",
                        help="relative_to_max: TPM / gene's top transcript TPM per tissue; "
                             "absolute: mean TPM. Default: %(default)s")
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    table = tissue_mean_tpm(args.input)
    if args.scaling == "relative_to_max":
        table = relative_to_max(table)
    table.columns = list(table.columns[:2]) + [tissue + "_TPM" for tissue in table.columns[2:]]
    table.to_csv(args.output, sep="\t", index=False)
    print("Wrote %d transcripts x %d tissues (%s) -> %s"
          % (len(table), len(table.columns) - 2, args.scaling, args.output))