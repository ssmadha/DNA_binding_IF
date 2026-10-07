#!/usr/bin/env python
"""
Fit the ASIF impact-factor sigmoid to experimental isoform interaction data.

compute_asif.py scores a transcript whose domains have coverages c_1..c_n as

    impact_factor = 1 - mean_i( kept(c_i) )
    kept(c) = 1 if c == 1 (fully intact), else sigmoid(alpha * (c - beta))

so 1 - impact_factor is the predicted fraction of function kept. This script
fits alpha and beta (one set per assay/domain type) to Lambourne et al.
2025, Mol Cell 85:1445 (PMID
40147441), supplementary tables:

    Table_S1.tsv  - isoform clones (isoform_status marks each gene's
                    reference; ensembl_transcript_ids, aa_seq)
    Table_S3.tsv  - eY1H protein-DNA interactions (PDIs) per clone and DNA
                    bait (TRUE = PDI, FALSE = verified non-PDI, blank = not
                    successfully tested)
    Table_S5.tsv  - Y2H protein-protein interactions (PPIs) per clone and
                    partner (Y2H_result TRUE/FALSE/NA)
    Table_S13.tsv - DBD_pct_lost_in_alt and Ensembl gene ID per pair

Each alternative isoform is compared with its gene's reference on the baits
or partners tested for both; y = fraction of the reference's interactions the
alternative keeps (pairs where the reference has >= 1).

--assay pdi: one domain per pair, the DBD; coverage = 1 - DBD_pct_lost_in_alt / 100.
--assay ppi: the reference's PPI domains as the pipeline classifies them
    (InterPro domains of type PPI, from the gene's Ensembl protein closest to
    the reference clone, carried onto the clone by alignment; overlapping
    domains merged as in matching._merge_overlapping, including DNA-binding +
    PPI domains suppressing overlapping PPI domains). Coverage of each = the
    fraction of its residues aligned to an identical residue in the
    alternative's sequence.

alpha and beta are fitted by weighted least squares (each pair equally, or with --weighting
interactions by its number of reference interactions) over an alpha x beta
grid. Near-step data leave the fit almost flat in alpha, so besides the
best fit the script reports the least steep alpha whose SSE is within
SSE_TOLERANCE of it ("chosen"), plus an alpha profile; pairs are bootstrapped
for confidence intervals of the chosen fit.

Outputs <prefix>_pairs.tsv (per-pair counts, category - matching Table
S13's PDI/PPI_category - and domain coverages) and <prefix>.png.
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Align
from scipy.optimize import minimize_scalar

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "bin"))
from compute_asif import domain_retained  # noqa: E402

# The data can be close to a step, so alpha is capped; beta is in [0, 1].
ALPHA_MAX = 200
ALPHA_GRID = np.unique(np.round(np.geomspace(0.5, ALPHA_MAX, 120), 2))
BETA_GRID = np.linspace(0, 1, 501)
# "Chosen" fit: the smallest alpha whose SSE is within this fraction of the best.
SSE_TOLERANCE = 0.01
PROFILE_ALPHAS = (5, 10, 20, 30, 40, 60, 80, 120, 200)
COMPARISON_CURVES = {"notebook (6, 0.5)": (6, 0.5), "pipeline default (63, 0.3)": (63, 0.3)}
ASSAY_LABELS = {"pdi": ("PDI", "eY1H baits", "DBD"), "ppi": ("PPI", "Y2H partners", "PPI domains")}

# Reference palette: categorical slots 1-3, neutral ink for data points.
FIT_COLOR = "#2a78d6"
COMPARISON_COLORS = ["#eb6834", "#1baf7a"]
POINT_COLOR = "#52514e"
TEXT_COLOR = "#0b0b0b"
SURFACE = "#fcfcfb"


class PairCoverages:
    """Per-pair domain coverages, flattened for vectorized prediction."""

    def __init__(self, coverages):
        self.coverages = [np.asarray(values, dtype=float) for values in coverages]
        self.counts = np.array([len(values) for values in self.coverages])
        self.values = np.concatenate(self.coverages)
        self.starts = np.concatenate([[0], np.cumsum(self.counts)[:-1]])

    def subset(self, index):
        return PairCoverages([self.coverages[i] for i in index])

    def predict(self, alpha, beta):
        """
        Per pair, mean domain_retained over its domains (= 1 - impact
        factor). beta may be an array; returns one row per beta.
        """
        betas = np.atleast_1d(beta)[:, None]
        kept = domain_retained(self.values[None, :], alpha, betas)
        return np.add.reduceat(kept, self.starts, axis=1) / self.counts


def category(pair, label):
    """Change in interactions from reference to alternative, as in Table S13."""
    if pair.n_tested_both == 0:
        return "not tested in both"
    if pair.ref_interactions == 0 and pair.alt_interactions == 0:
        return "neither interacts"
    if pair.alt_interactions == 0:
        return label + " loss"
    if pair.ref_interactions == 0:
        return label + " gain"
    if pair.ref_only == 0 and pair.alt_only == 0:
        return "no " + label + " change"
    return label + " rewire"


def compare_isoforms(calls, s1, label):
    """
    Compare each alternative isoform with its gene's reference.

    Parameters
    ----------
    calls: pandas.DataFrame
        Clone ID index, one column per bait/partner; True/False/NaN.
    s1: pandas.DataFrame
        Table S1.
    label: str
        "PDI" or "PPI", for the category names.

    Returns
    -------
    pandas.DataFrame
        One row per pair with both isoforms in calls.
    """
    gene_of = s1.set_index("clone_id")["gene_symbol"]
    status = s1.set_index("clone_id")["isoform_status"]
    references = s1[s1["isoform_status"].str.contains("reference")].set_index("gene_symbol")["clone_id"]
    rows = []
    for alternative in calls.index:
        gene = gene_of.get(alternative)
        reference = references.get(gene)
        if "alternative" not in str(status.get(alternative)) or reference not in calls.index:
            continue
        tested = calls.loc[reference].notna() & calls.loc[alternative].notna()
        ref_hits = calls.loc[reference][tested].astype(bool)
        alt_hits = calls.loc[alternative][tested].astype(bool)
        rows.append({"gene_symbol": gene, "reference_isoform": reference,
                     "alternative_isoform": alternative,
                     "n_tested_both": int(tested.sum()),
                     "ref_interactions": int(ref_hits.sum()), "alt_interactions": int(alt_hits.sum()),
                     "shared": int((ref_hits & alt_hits).sum()),
                     "ref_only": int((ref_hits & ~alt_hits).sum()),
                     "alt_only": int((~ref_hits & alt_hits).sum())})
    pairs = pd.DataFrame(rows)
    pairs["category"] = pairs.apply(category, axis=1, label=label)
    pairs["frac_ref_retained"] = pairs["shared"] / pairs["ref_interactions"].where(pairs["ref_interactions"] > 0)
    return pairs


def load_calls(supp_dir, assay):
    """Interaction calls (clone x bait/partner) for the assay."""
    if assay == "pdi":
        s3 = pd.read_csv(os.path.join(supp_dir, "Table_S3.tsv"), sep="\t").set_index("clone_id")
        return s3.drop(columns="gene_symbol")
    s5 = pd.read_csv(os.path.join(supp_dir, "Table_S5.tsv"), sep="\t").dropna(subset=["Y2H_result"])
    s5["Y2H_result"] = s5["Y2H_result"].astype(str).str.lower() == "true"
    return s5.pivot_table(index="ad_clone_id", columns="db_gene_symbol", values="Y2H_result", aggfunc="max")


def make_aligner():
    """Global protein aligner with cheap gap extension and free end gaps, for isoforms."""
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 5
    aligner.mismatch_score = -4
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.end_gap_score = 0
    return aligner


def identical_residue_map(aligner, query, target):
    """
    Map query residue positions to target positions where the aligned
    residues are identical.

    Returns
    -------
    dict[int, int]
    """
    alignment = aligner.align(query, target)[0]
    mapping = {}
    for (q_start, q_end), (t_start, t_end) in zip(*alignment.aligned):
        for offset in range(q_end - q_start):
            if query[q_start + offset] == target[t_start + offset]:
                mapping[q_start + offset] = t_start + offset
    return mapping


def ppi_domain_coverages(pairs, s1, gene_ids, reference_dir):
    """
    Coverage of each reference's PPI domains by its alternative.

    Parameters
    ----------
    pairs: pandas.DataFrame
        Output of compare_isoforms.
    s1: pandas.DataFrame
        Table S1.
    gene_ids: dict
        Gene symbol -> Ensembl gene ID.
    reference_dir: str
        Directory holding the pipeline's reference files.

    Returns
    -------
    pandas.DataFrame
        Per alternative: source_protein, source_identity (fraction of the
        reference clone's residues identical to the source Ensembl
        protein), ppi_domains (";"-separated IDs) and coverages (list).
    """
    from TF_ASIF.gene import Gene
    from TF_ASIF.matching import _merge_overlapping
    from TF_ASIF.transcript import Transcript

    files = {name: os.path.join(reference_dir, path) for name, path in {
        "gtf": "Homo_sapiens.GRCh38.109.gtf.gz", "cds": "Homo_sapiens.GRCh38.cds.all.fa.gz",
        "interpro": "Homo_sapiens.GRCh38.interpro_domains.tsv.gz",
        "pfam": "Homo_sapiens.GRCh38.pfam_domains.tsv.gz",
        "entry_types": "interpro_90.0_entry.list"}.items()}
    gtf_index = Gene._get_gtf_index(files["gtf"])
    cds = Transcript._get_cds_sequences(files["cds"])
    clones = s1.set_index("clone_id")
    aligner = make_aligner()

    reference_domains = {}
    for reference in pairs["reference_isoform"].unique():
        clone = clones.loc[reference]
        ref_seq = clone["aa_seq"].rstrip("*")
        listed = str(clone["ensembl_transcript_ids"]).split("|")
        candidates = []
        for isoform in gtf_index.get(gene_ids.get(clone["gene_symbol"]), []):
            if isoform["protein_id"] is None or isoform["id"] not in cds or len(cds[isoform["id"]]) % 3:
                continue
            protein = str(cds[isoform["id"]].translate()).rstrip("*")
            mapping = identical_residue_map(aligner, protein, ref_seq)
            candidates.append((len(mapping) / len(ref_seq), isoform["id"] in listed, isoform["protein_id"], mapping))
        if not candidates:
            reference_domains[reference] = (None, np.nan, [])
            continue
        identity, _, ensp, mapping = max(candidates, key=lambda candidate: candidate[:2])
        mapped = []
        for domain in Transcript.interpro_domain_objects(ensp, files["interpro"], files["pfam"], files["entry_types"]):
            if not {"DNA-binding", "PPI"} & set(domain.types):
                continue
            residues = sorted({mapping[residue] for part in domain.pos.parts
                               for residue in range(part.start, part.end) if residue in mapping})
            if residues:
                mapped.append((domain, ensp, [(residue,) for residue in residues]))
        kept = [(entry[0].domain_id, {codon[0] for codon in entry[2]})
                for entry in _merge_overlapping(mapped) if "PPI" in entry[0].types]
        reference_domains[reference] = (ensp, identity, kept)

    rows = []
    for _, pair in pairs.iterrows():
        ensp, identity, domains = reference_domains[pair["reference_isoform"]]
        alt_seq = clones.loc[pair["alternative_isoform"], "aa_seq"].rstrip("*")
        ref_seq = clones.loc[pair["reference_isoform"], "aa_seq"].rstrip("*")
        kept_residues = set(identical_residue_map(aligner, ref_seq, alt_seq))
        rows.append({"alternative_isoform": pair["alternative_isoform"], "source_protein": ensp,
                     "source_identity": identity,
                     "ppi_domains": ";".join(domain_id for domain_id, _ in domains),
                     "coverages": [len(residues & kept_residues) / len(residues) for _, residues in domains]})
    return pd.DataFrame(rows)


def load_pairs(supp_dir, assay, reference_dir):
    """
    Reference/alternative pairs with interaction counts and per-domain coverages.

    Returns
    -------
    pandas.DataFrame
        compare_isoforms columns plus DBD_pct_lost_in_alt, coverages (list
        of per-domain coverages, empty if none) and mean_coverage.
    """
    s1 = pd.read_csv(os.path.join(supp_dir, "Table_S1.tsv"), sep="\t")
    s13 = pd.read_csv(os.path.join(supp_dir, "Table_S13.tsv"), sep="\t")
    pairs = compare_isoforms(load_calls(supp_dir, assay), s1, ASSAY_LABELS[assay][0])
    pairs = pairs.merge(s13[["reference_isoform", "alternative_isoform", "DBD_pct_lost_in_alt"]],
                        on=["reference_isoform", "alternative_isoform"], how="left")
    if assay == "pdi":
        pairs["coverages"] = [[] if pd.isna(lost) else [1 - lost / 100] for lost in pairs["DBD_pct_lost_in_alt"]]
    else:
        gene_ids = dict(zip(s13["gene_symbol"], s13["Ensembl_gene_ID"]))
        pairs = pairs.merge(ppi_domain_coverages(pairs, s1, gene_ids, reference_dir), on="alternative_isoform")
    pairs["n_domains"] = pairs["coverages"].map(len)
    pairs["mean_coverage"] = pairs["coverages"].map(lambda values: np.mean(values) if values else np.nan)
    return pairs


def alpha_profile(coverages, retained, weights, alphas=ALPHA_GRID):
    """
    Best beta and weighted SSE at each fixed alpha (beta grid search, then
    refined). A flat SSE means the data don't determine alpha.

    Returns
    -------
    numpy.ndarray
        Rows of (alpha, beta, sse), in the order of alphas.
    """
    def sse(alpha, beta):
        return np.sum(weights * (retained - coverages.predict(alpha, beta)) ** 2, axis=1)

    rows = []
    for alpha in alphas:
        grid_sse = sse(alpha, BETA_GRID)
        i = int(grid_sse.argmin())
        low, high = BETA_GRID[max(i - 1, 0)], BETA_GRID[min(i + 1, len(BETA_GRID) - 1)]
        result = minimize_scalar(lambda beta: sse(alpha, beta)[0], bounds=(low, high), method="bounded")
        rows.append((alpha, result.x, result.fun) if result.fun < grid_sse[i] else (alpha, BETA_GRID[i], grid_sse[i]))
    return np.array(rows)


def fit_sigmoid(coverages, retained, weights):
    """
    Weighted least-squares fit of alpha and beta.

    Returns
    -------
    tuple
        (best (alpha, beta, sse), chosen (alpha, beta, sse) = the smallest
        alpha within SSE_TOLERANCE of the best SSE)
    """
    profile = alpha_profile(coverages, retained, weights)
    best = profile[profile[:, 2].argmin()]
    chosen = profile[profile[:, 2] <= best[2] * (1 + SSE_TOLERANCE) + 1e-12][0]
    return tuple(best), tuple(chosen)


def bootstrap(coverages, retained, weights, n_boot, seed):
    """
    Refit the chosen alpha/beta on n_boot resamples of the pairs. Resamples
    where no pair keeps any interaction can't place the curve and are
    skipped.

    Parameters
    ----------
    coverages: PairCoverages
    retained, weights: numpy.ndarray
        Per-pair fraction of reference interactions kept, and weight.
    n_boot: int
    seed: int

    Returns
    -------
    tuple
        (numpy.ndarray of chosen (alpha, beta) fits, number of skipped resamples)
    """
    rng = np.random.default_rng(seed)
    fits = []
    skipped = 0
    for _ in range(n_boot):
        index = rng.integers(0, len(retained), len(retained))
        if not retained[index].any():
            skipped += 1
            continue
        _, chosen = fit_sigmoid(coverages.subset(index), retained[index], weights[index])
        fits.append(chosen[:2])
    return np.array(fits), skipped


def plot_fit(fitted, alpha, beta, assay, output_png):
    """Data points, fitted curve and comparison curves; full range and a zoom near full coverage."""
    _, interactions, domains = ASSAY_LABELS[assay]
    x = np.linspace(0, 1, 1001)
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.5), facecolor=SURFACE,
                                gridspec_kw={"width_ratios": [1.3, 1]})
    for axis, x_range in zip(axes, [(-0.02, 1.02), (0.75, 1.005)]):
        axis.set_facecolor(SURFACE)
        sizes = 12 + 4 * fitted["ref_interactions"].clip(upper=40)
        axis.scatter(fitted["mean_coverage"], fitted["frac_ref_retained"], s=sizes,
                     color=POINT_COLOR, alpha=0.45, edgecolors=SURFACE, linewidths=1,
                     label="isoform pair (size ∝ ref. %s)" % interactions.split()[1], zorder=3)
        axis.plot(x, domain_retained(x, alpha, beta), color=FIT_COLOR, lw=2,
                  label="fit (%.1f, %.3f)" % (alpha, beta), zorder=4)
        for (name, (a, b)), color in zip(COMPARISON_CURVES.items(), COMPARISON_COLORS):
            axis.plot(x, domain_retained(x, a, b), color=color, lw=2, ls="--", label=name, zorder=2)
        axis.set_xlim(*x_range)
        axis.set_ylim(-0.04, 1.04)
        axis.grid(color="#e4e3df", lw=0.8)
        axis.set_axisbelow(True)
        for side in ("top", "right"):
            axis.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            axis.spines[side].set_color("#b5b4ae")
        axis.tick_params(colors="#52514e")
        axis.set_xlabel("%scoverage of reference %s in alternative"
                        % ("" if assay == "pdi" else "Mean ", domains), color=TEXT_COLOR)
    axes[0].set_ylabel("Fraction of reference %s\nalso bound by alternative" % interactions, color=TEXT_COLOR)
    axes[0].set_title("All pairs (n = %d)" % len(fitted), color=TEXT_COLOR, loc="left")
    axes[1].set_title("Zoom: coverage 0.75-1", color=TEXT_COLOR, loc="left")
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(handles, labels, frameon=False, loc="lower center", ncol=4, fontsize=9, labelcolor=TEXT_COLOR)
    figure.tight_layout(rect=(0, 0.07, 1, 1))
    figure.savefig(output_png, dpi=200, facecolor=SURFACE)


def get_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-a", "--assay", choices=["pdi", "ppi"], default="pdi",
                        help="pdi: eY1H vs DBD coverage; ppi: Y2H vs PPI-domain coverage. Default: %(default)s")
    parser.add_argument("-d", "--supp-dir", default="Gloria_paper_supplemental",
                        help="Directory with Table_S1/S3/S5/S13.tsv. Default: %(default)s")
    parser.add_argument("-r", "--reference-dir", default="reference_data",
                        help="Pipeline reference files (GTF, CDS FASTA, InterPro/Pfam domains, "
                             "InterPro entry types), for --assay ppi. Default: %(default)s")
    parser.add_argument("-o", "--output-prefix", required=True,
                        help="Writes <prefix>_pairs.tsv and <prefix>.png")
    parser.add_argument("-w", "--weighting", choices=["pairs", "interactions"], default="pairs",
                        help="Weight each pair equally, or by its number of reference interactions. "
                             "Default: %(default)s")
    parser.add_argument("-x", "--exclude-full-coverage", action="store_true",
                        help="Fit alpha/beta only on pairs whose alternative lost part of a domain "
                             "(mean coverage < 1)")
    parser.add_argument("-n", "--n-bootstrap", type=int, default=1000,
                        help="Bootstrap resamples for confidence intervals. Default: %(default)s")
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    pairs = load_pairs(args.supp_dir, args.assay, args.reference_dir)
    pairs.to_csv(args.output_prefix + "_pairs.tsv", sep="\t", index=False)

    def pair_weights(subset):
        return (subset["ref_interactions"] if args.weighting == "interactions"
                else pd.Series(1, index=subset.index)).to_numpy(float)

    with_interactions = pairs[pairs["ref_interactions"] > 0]
    eligible = with_interactions[with_interactions["n_domains"] > 0]
    fitted =eligible[eligible["mean_coverage"] < 1] if args.exclude_full_coverage else eligible
    coverages = PairCoverages(fitted["coverages"])
    retained = fitted["frac_ref_retained"].to_numpy()
    weights = pair_weights(fitted)

    best, chosen = fit_sigmoid(coverages, retained, weights)
    alpha, beta, _ = chosen
    boot, skipped = bootstrap(coverages, retained, weights, args.n_bootstrap, args.seed)
    alpha_ci = np.percentile(boot[:, 0], [2.5, 97.5])
    beta_ci = np.percentile(boot[:, 1], [2.5, 97.5])

    print("Pairs compared: %d; reference has >= 1 interaction: %d; with >= 1 domain: %d; fitted: %d from %d genes"
          % (len(pairs), len(with_interactions), len(eligible), len(fitted), fitted["gene_symbol"].nunique()))
    print("Categories: %s" % pairs["category"].value_counts().to_dict())
    print("Domains per fitted pair: %s" % fitted["n_domains"].value_counts().sort_index().to_dict())
    print("Weighting: %s" % args.weighting)
    print("Best fit:   alpha = %.2f, beta = %.3f, SSE %.4f" % best)
    print("Chosen (smallest alpha within %g%% of best SSE): alpha = %.2f (95%% CI %.2f-%.2f), "
          "beta = %.3f (95%% CI %.3f-%.3f), SSE %.4f"
          % (100 * SSE_TOLERANCE, alpha, *alpha_ci, beta, *beta_ci, chosen[2]))
    print("Bootstrap: %d fits, %d resamples skipped (no interaction kept); %d at alpha bound %d"
          % (len(boot), skipped, (boot[:, 0] > ALPHA_MAX - 1).sum(), ALPHA_MAX))
    residual = retained - coverages.predict(alpha, beta)[0]
    print("Weighted RMSE: chosen fit %.3f" % np.sqrt(np.average(residual ** 2, weights=weights)), end="")
    for name, (a, b) in COMPARISON_CURVES.items():
        residual = retained - coverages.predict(a, b)[0]
        print("; %s %.3f" % (name, np.sqrt(np.average(residual ** 2, weights=weights))), end="")
    print()
    print("Alpha profile (if SSE is flat over a range, alpha is not determined by the data):")
    print("  alpha    beta     SSE  kept@0.99  @0.95  @0.90  @0.50")
    for profile_alpha, profile_beta, sse in alpha_profile(coverages, retained, weights, PROFILE_ALPHAS):
        print("  %5g  %.4f  %6.4f  %9.2f  %5.2f  %5.2f  %5.2f"
              % (profile_alpha, profile_beta, sse,
                 *domain_retained(np.array([0.99, 0.95, 0.9, 0.5]), profile_alpha, profile_beta)))
    plot_fit(fitted, alpha, beta, args.assay, args.output_prefix + ".png")
    print("Wrote %s_pairs.tsv and %s.png" % (args.output_prefix, args.output_prefix))