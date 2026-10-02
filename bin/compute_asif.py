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
results: merged results TSV from merge_gene_results.py, either
    - one row per transcript and domain with a "coverage" column (current
      download_gene.py output), or
    - one row per transcript with n_domains and a "[a, b, ...]"
      domain_coverage list (older runs).
expression: TSV with gene_id, transcript_id, then one column per
    tissue/sample named "<tissue><suffix>" (suffix "_TPM" by default).
    Columns not ending in the suffix are ignored.

Output columns: gene_id, transcript_id, n_domains, domain_coverage, then
(when the results have a domain_type column) n_<type>_domains and
<type>_coverage for each domain type in DOMAIN_TYPES, then
"<tissue>_tpm", "<tissue>_asif" for each tissue, in the expression file's
column order. A domain of both types counts under each. Only transcripts
present in both inputs are kept.
"""
import argparse
import sys

import numpy as np
import pandas as pd

KEY_COLUMNS = ["gene_id", "transcript_id"]
# domain_type value (";"-separated in the results) -> output column prefix.
DOMAIN_TYPES = {"DNA-binding": "dna_binding", "PPI": "ppi"}


def coverage_list(values):
    """Format coverages as a "[a, b, ...]" list string."""
    return str([float(value) for value in values])


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


def to_per_transcript(results_df):
    """
    Collapse a one-row-per-(transcript, domain) results table into one
    row per transcript with n_domains and a "[a, b, ...]"
    domain_coverage list, plus, if the table has a domain_type column,
    n_<type>_domains and <type>_coverage for each DOMAIN_TYPES type
    (a domain of both types counts under each). Tables already in that
    form are returned unchanged. Domains with no mappable residues (NaN
    coverage) are left out.

    Parameters
    ----------
    results_df: pandas.DataFrame

    Returns
    -------
    pandas.DataFrame
    """
    if "domain_coverage" in results_df.columns:
        return results_df
    scored = results_df.dropna(subset=["coverage"])
    grouped = scored.groupby(KEY_COLUMNS, sort=False)["coverage"]
    per_transcript = grouped.agg(n_domains="size", domain_coverage=coverage_list).reset_index()
    if "domain_type" not in scored.columns:
        return per_transcript
    types = scored["domain_type"].fillna("").str.split(";")
    for domain_type, prefix in DOMAIN_TYPES.items():
        of_type = scored[types.map(lambda value: domain_type in value)]
        by_type = of_type.groupby(KEY_COLUMNS, sort=False)["coverage"].agg(
            **{"n_%s_domains" % prefix: "size", "%s_coverage" % prefix: coverage_list}).reset_index()
        per_transcript = per_transcript.merge(by_type, on=KEY_COLUMNS, how="left")
        per_transcript["n_%s_domains" % prefix] = per_transcript["n_%s_domains" % prefix].fillna(0).astype(int)
        per_transcript["%s_coverage" % prefix] = per_transcript["%s_coverage" % prefix].fillna("[]")
    return per_transcript


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

    results_df = to_per_transcript(results_df)
    merged = results_df.merge(expression_df[KEY_COLUMNS + expression_columns], on=KEY_COLUMNS)
    factors = merged["domain_coverage"].map(lambda value: impact_factor(parse_coverage(value), alpha, beta))

    type_columns = [column for prefix in DOMAIN_TYPES.values()
                    for column in ("n_%s_domains" % prefix, "%s_coverage" % prefix) if column in merged.columns]
    out = merged[["gene_id", "transcript_id", "n_domains", "domain_coverage"] + type_columns].copy()
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

    n_results = len(to_per_transcript(results_df))
    if asif_df.empty:
        print("WARNING: no transcripts in %s matched %s" % (args.results, args.expression), file=sys.stderr)
    print("Wrote %d of %d transcripts (%d not in expression file) -> %s"
          % (len(asif_df), n_results, n_results - len(asif_df), args.output), file=sys.stderr)