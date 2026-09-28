#!/usr/bin/env python
"""
Compute per-tissue alternative-splicing impact factors (ASIF) from merged
domain-coverage results and a transcript expression table.

For each transcript with per-domain coverages c_1..c_n (from
merge_gene_results.py):

    impact_factor = 1 - mean_i( sigmoid(alpha * (c_i - beta)) )

and for each tissue/sample column:

    <tissue>_asif = impact_factor * <tissue>_tpm

Inputs
------
results: TSV with gene_id, transcript_id, n_domains, domain_coverage
    (domain_coverage as a "[a, b, ...]" list).
expression: TSV with gene_id, transcript_id, then one column per
    tissue/sample named "<tissue><suffix>" (suffix "_TPM" by default).
    Columns not ending in the suffix are ignored.

Output columns: gene_id, transcript_id, n_domains, domain_coverage, then
"<tissue>_tpm", "<tissue>_asif" for each tissue, in the expression file's
column order. Only transcripts present in both inputs are kept.
"""
import argparse
import sys

import numpy as np
import pandas as pd

KEY_COLUMNS = ["gene_id", "transcript_id"]


def parse_coverage(coverage_str):
    """
    Parse a "[a, b, ...]" coverage list string.

    Parameters
    ----------
    coverage_str: str

    Returns
    -------
    numpy.ndarray
    """
    inner = coverage_str.strip()[1:-1].strip()
    if not inner:
        return np.array([], dtype=float)
    return np.array([float(value) for value in inner.split(",")])


def impact_factor(coverages, alpha, beta):
    """
    Impact factor of one transcript.

    Parameters
    ----------
    coverages: numpy.ndarray
        Per-domain coverage fractions of this transcript.
    alpha: float
        Sigmoid steepness.
    beta: float
        Coverage at which the sigmoid is centred.

    Returns
    -------
    float
        1 - mean sigmoid(alpha * (coverage - beta)), or NaN if there are
        no domains.
    """
    if coverages.size == 0:
        return np.nan
    return 1 - np.mean(1 / (1 + np.exp(-alpha * (coverages - beta))))


def compute_asif(results_df, expression_df, alpha, beta, suffix="_TPM"):
    """
    Build the ASIF table.

    Parameters
    ----------
    results_df: pandas.DataFrame
        Merged coverage results.
    expression_df: pandas.DataFrame
        Transcript expression table.
    alpha: float
    beta: float
    suffix: str
        Suffix identifying expression columns in expression_df.

    Returns
    -------
    pandas.DataFrame
    """
    expression_columns = [column for column in expression_df.columns
                          if column not in KEY_COLUMNS and column.endswith(suffix)]
    if not expression_columns:
        raise ValueError("No expression columns ending in '%s' found" % suffix)

    merged = results_df.merge(expression_df[KEY_COLUMNS + expression_columns], on=KEY_COLUMNS)
    factors = merged["domain_coverage"].map(lambda value: impact_factor(parse_coverage(value), alpha, beta))

    out = merged[["gene_id", "transcript_id", "n_domains", "domain_coverage"]].copy()
    new_columns = {}
    for column in expression_columns:
        tissue = column[:-len(suffix)]
        new_columns[tissue + "_tpm"] = merged[column]
        new_columns[tissue + "_asif"] = factors * merged[column]
    return pd.concat([out, pd.DataFrame(new_columns, index=merged.index)], axis=1)


def get_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-r", "--results", required=True,
                        help="Merged results TSV from merge_gene_results.py")
    parser.add_argument("-e", "--expression", required=True,
                        help="Transcript expression TSV (gene_id, transcript_id, <tissue>_TPM, ...)")
    parser.add_argument("-o", "--output", required=True,
                        help="Path to write the ASIF table to")
    parser.add_argument("-a", "--alpha", type=float, default=63,
                        help="Sigmoid steepness. Default: %(default)s")
    parser.add_argument("-b", "--beta", type=float, default=0.3,
                        help="Coverage at which the sigmoid is centred. Default: %(default)s")
    parser.add_argument("-s", "--suffix", default="_TPM",
                        help="Suffix of expression columns in the expression file. Default: %(default)s")
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    results_df = pd.read_csv(args.results, sep="\t")
    expression_df = pd.read_csv(args.expression, sep="\t")
    asif_df = compute_asif(results_df, expression_df, args.alpha, args.beta, args.suffix)
    asif_df.to_csv(args.output, sep="\t", index=False)

    n_results = len(results_df)
    if asif_df.empty:
        print("WARNING: no transcripts in %s matched %s" % (args.results, args.expression), file=sys.stderr)
    print("Wrote %d of %d transcripts (%d not in expression file) -> %s"
          % (len(asif_df), n_results, n_results - len(asif_df), args.output), file=sys.stderr)