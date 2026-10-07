#!/usr/bin/env python
"""
Compute per-tissue alternative-splicing impact factors (ASIF) from merged
domain-coverage results and a transcript expression table.

For each transcript with per-domain coverages c_1..c_n (from
merge_gene_results.py), each domain keeps a fraction of its function

    kept(c) = 1                            if c == 1 (fully intact)
            = sigmoid(alpha * (c - beta))  otherwise

with alpha and beta set separately for DNA-binding and PPI domains. A domain of both
types, or of no recorded type (alignment mode, older runs), uses the
parameters chosen by --both-type (DNA-binding by default). Then

    impact_factor = 1 - mean_i( kept(c_i) )

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
# domain_type value (";"-separated in the results) -> output column prefix,
# also the name of that type's parameter set.
DOMAIN_TYPES = {"DNA-binding": "dna_binding", "PPI": "ppi"}
# (alpha, beta) per type, from code/scripts/fit_asif_sigmoid.py -x on the
# partially covered domains of Lambourne et al. 2025: PPI (Y2H) = best fit;
# DNA-binding (eY1H) = least steep alpha within 1% of the best fit, since the
# fit is flat for alpha >= ~40.
DEFAULT_PARAMS = {"dna_binding": (48.8, 0.958), "ppi": (17.0, 0.965)}


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


def sigmoid(x):
    """Logistic function, without overflow for large |x|."""
    return np.exp(-np.logaddexp(0, -x))


def domain_retained(coverages, alpha, beta):
    """
    Fraction of a domain's function kept at the given coverage(s).

    Parameters
    ----------
    coverages: numpy.ndarray or float
        Domain coverage fractions.
    alpha: float
        Sigmoid steepness.
    beta: float
        Coverage at which the sigmoid is centred.

    Returns
    -------
    numpy.ndarray or float
        1 for a fully covered domain, else sigmoid(alpha * (coverage - beta)).
    """
    coverages = np.asarray(coverages, dtype=float)
    return np.where(coverages >= 1, 1.0, sigmoid(alpha * (coverages - beta)))


def impact_factor(coverages, alpha, beta):
    """
    Impact factor of one transcript whose domains share one parameter set.

    Parameters
    ----------
    coverages: numpy.ndarray
        Per-domain coverage fractions of this transcript.
    alpha, beta: float
        See domain_retained.

    Returns
    -------
    float
        1 - mean domain_retained, or NaN if there are no domains.
    """
    if coverages.size == 0:
        return np.nan
    return 1 - np.mean(domain_retained(coverages, alpha, beta))


def parameter_type(domain_type, both_type):
    """
    Parameter set ("dna_binding" or "ppi") for a domain_type value
    (";"-separated); both types or none -> both_type.
    """
    types = {DOMAIN_TYPES[value] for value in str(domain_type).split(";") if value in DOMAIN_TYPES}
    return types.pop() if len(types) == 1 else both_type


def transcript_impact_factors(results_df, params, both_type):
    """
    Impact factor of every transcript in the results.

    Parameters
    ----------
    results_df: pandas.DataFrame
        Merged coverage results (one row per transcript and domain, or the
        older one row per transcript with a domain_coverage list, whose
        domains are untyped).
    params: dict
        Parameter type -> (alpha, beta).
    both_type: str
        Parameter type for domains of both types or no recorded type.

    Returns
    -------
    pandas.DataFrame
        gene_id, transcript_id, impact_factor.
    """
    if "domain_coverage" in results_df.columns:
        factors = results_df[KEY_COLUMNS].copy()
        factors["impact_factor"] = results_df["domain_coverage"].map(
            lambda value: impact_factor(parse_coverage(value), *params[both_type]))
        return factors
    scored = results_df.dropna(subset=["coverage"])
    domain_types = scored["domain_type"].fillna("") if "domain_type" in scored.columns else pd.Series("", index=scored.index)
    types = domain_types.map(lambda value: parameter_type(value, both_type))
    kept = pd.Series(np.nan, index=scored.index)
    for param_type, (alpha, beta) in params.items():
        of_type = types == param_type
        kept[of_type] = domain_retained(scored.loc[of_type, "coverage"].to_numpy(float), alpha, beta)
    kept_mean = kept.groupby([scored["gene_id"], scored["transcript_id"]], sort=False).mean()
    return (1 - kept_mean).rename("impact_factor").reset_index()


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


def compute_asif(results_df, expression_df, params=None, both_type="dna_binding", suffix="_TPM"):
    """
    Build the ASIF table.

    Parameters
    ----------
    results_df: pandas.DataFrame
        Merged coverage results.
    expression_df: pandas.DataFrame
        Transcript expression table.
    params: dict, optional
        Parameter type ("dna_binding", "ppi") -> (alpha, beta).
        Defaults to DEFAULT_PARAMS.
    both_type: str
        Parameter type for domains of both types or no recorded type.
    suffix: str
        Suffix identifying expression columns in expression_df.

    Returns
    -------
    pandas.DataFrame
    """
    if params is None:
        params = DEFAULT_PARAMS
    expression_columns = [column for column in expression_df.columns
                          if column not in KEY_COLUMNS and column.endswith(suffix)]
    if not expression_columns:
        raise ValueError("No expression columns ending in '%s' found" % suffix)

    per_transcript = to_per_transcript(results_df).merge(
        transcript_impact_factors(results_df, params, both_type), on=KEY_COLUMNS)
    merged = per_transcript.merge(expression_df[KEY_COLUMNS + expression_columns], on=KEY_COLUMNS)
    factors = merged["impact_factor"]

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
    for param_type, (alpha, beta) in DEFAULT_PARAMS.items():
        option = param_type.replace("_", "-")
        parser.add_argument("--%s-alpha" % option, type=float, default=alpha,
                            help="Sigmoid steepness for %s domains. Default: %%(default)s" % param_type)
        parser.add_argument("--%s-beta" % option, type=float, default=beta,
                            help="Coverage at which the %s sigmoid is centred. Default: %%(default)s" % param_type)
    parser.add_argument("-t", "--both-type", choices=list(DEFAULT_PARAMS), default="dna_binding",
                        help="Parameters for domains of both types or of no recorded type. Default: %(default)s")
    parser.add_argument("-s", "--suffix", default="_TPM",
                        help="Suffix of expression columns in the expression file. Default: %(default)s")
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    results_df = pd.read_csv(args.results, sep="\t")
    expression_df = pd.read_csv(args.expression, sep="\t")
    params = {param_type: (getattr(args, param_type + "_alpha"), getattr(args, param_type + "_beta"))
              for param_type in DEFAULT_PARAMS}
    asif_df = compute_asif(results_df, expression_df, params, args.both_type, args.suffix)
    asif_df.to_csv(args.output, sep="\t", index=False)

    n_results = len(to_per_transcript(results_df))
    if asif_df.empty:
        print("WARNING: no transcripts in %s matched %s" % (args.results, args.expression), file=sys.stderr)
    print("Wrote %d of %d transcripts (%d not in expression file) -> %s"
          % (len(asif_df), n_results, n_results - len(asif_df), args.output), file=sys.stderr)