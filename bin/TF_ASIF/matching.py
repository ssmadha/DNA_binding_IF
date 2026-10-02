"""
Segment-mode domain matching: place each domain / PPI binding site on the
transcript whose protein its positions refer to, map its residues to exact
codon bases in the gene's segment superisoform (see superisoform.py), and
score every transcript by the fraction of those residues it encodes.
"""
import warnings

import pandas as pd

from .domain import Domain
from .transcript import Transcript

# Columns of the per-transcript, per-domain coverage table written by
# download_gene.py (both match modes).
COVERAGE_COLUMNS = ["gene_id", "transcript_id", "domain_id", "domain_type", "source_transcript_id",
                    "positions", "n_residues", "n_covered", "coverage"]

_UNIPROT_DBS = ("Uniprot/SWISSPROT", "Uniprot/SPTREMBL", "Uniprot_isoform")
_uniprot_xrefs = {}


def load_uniprot_xrefs(uniprot_mapping_file):
    """
    Load and cache every UniProt xref of every transcript from an Ensembl
    protein-to-UniProt xref TSV file.

    Parameters
    ----------
    uniprot_mapping_file: str
        Path to a (optionally gzipped) Ensembl "*.uniprot.tsv.gz" file.

    Returns
    -------
    dict[str, list[tuple]]
        For each Ensembl Transcript ID, a list of (xref, db_name,
        identity) tuples, where identity is the lower of Ensembl's
        source/xref percent identities between the transcript's protein
        and the UniProt entry, or None if not given (Uniprot_isoform
        rows carry none).
    """
    if uniprot_mapping_file not in _uniprot_xrefs:
        mapping_df = pd.read_csv(uniprot_mapping_file, sep="\t", dtype=str)
        mapping_df = mapping_df[mapping_df["db_name"].isin(_UNIPROT_DBS)]
        xrefs = {}
        for transcript_id, xref, db_name, source_identity, xref_identity in zip(
                mapping_df["transcript_stable_id"], mapping_df["xref"], mapping_df["db_name"],
                mapping_df["source_identity"], mapping_df["xref_identity"]):
            try:
                identity = min(float(source_identity), float(xref_identity))
            except (TypeError, ValueError):
                identity = None
            xrefs.setdefault(transcript_id, []).append((xref, db_name, identity))
        _uniprot_xrefs[uniprot_mapping_file] = xrefs
    return _uniprot_xrefs[uniprot_mapping_file]


def choose_source_transcripts(gene_id, transcript_ids, uniprot_xrefs, accessions=None):
    """
    For each UniProt accession mapped to any of a gene's transcripts,
    pick the transcript whose protein best matches that accession's
    sequence: highest percent identity, then one tagged as isoform "-1",
    then the lowest transcript ID. Warns when the best match is below
    100% identity (a fallback).

    Parameters
    ----------
    gene_id: str
    transcript_ids: list[str]
    uniprot_xrefs: dict
        As returned by load_uniprot_xrefs.
    accessions: collection, optional
        Only consider these base accessions (e.g. those with PPI
        binding sites). Defaults to all.

    Returns
    -------
    dict[str, tuple[str, float | None]]
        (source transcript ID, identity) keyed by base UniProt
        accession (isoform suffix removed).
    """
    candidates = {}
    for transcript_id in sorted(transcript_ids):
        per_accession = {}
        for xref, db_name, identity in uniprot_xrefs.get(transcript_id, []):
            accession, _, isoform = xref.partition("-")
            if accessions is not None and accession not in accessions:
                continue
            best_identity, is_isoform_1 = per_accession.get(accession, (None, False))
            if identity is not None and (best_identity is None or identity > best_identity):
                best_identity = identity
            if db_name == "Uniprot_isoform" and isoform == "1":
                is_isoform_1 = True
            per_accession[accession] = (best_identity, is_isoform_1)
        for accession, (identity, is_isoform_1) in per_accession.items():
            score = (identity if identity is not None else -1, is_isoform_1)
            if accession not in candidates or score > candidates[accession][0]:
                candidates[accession] = (score, transcript_id, identity)

    sources = {}
    for accession, (_, transcript_id, identity) in candidates.items():
        if identity != 100:
            warnings.warn("%s: no transcript matches %s at 100%% identity; using %s (%s%%)"
                          % (gene_id, accession, transcript_id, identity))
        sources[accession] = (transcript_id, identity)
    return sources


def collect_domains(gene, superisoform, domain_filter, ppi_binding_site_file=None,
                    interpro_domains_file=None, uniprot_mapping_file=None, gtf_file=None,
                    pfam_domains_file=None, interpro_entry_types_file=None):
    """
    Gather a gene's classified domains, each paired with the transcript
    whose protein coordinates it is defined in: InterPro domains on
    their own Ensembl protein (for every transcript in the
    superisoform), PPI binding sites on their UniProt accession's
    source transcript (see choose_source_transcripts). Only InterPro
    domains of a requested classification are kept: DNA-binding for
    "dbi", PPI for "ppi_domain" (a domain of both types is kept under
    either, and reports both).

    Parameters
    ----------
    gene: Gene
    superisoform: Superisoform
    domain_filter: list
        Any of "ppi_domain", "dbi" (InterPro domains) and "ppi_bs"
        (PPI binding sites).
    ppi_binding_site_file, interpro_domains_file, uniprot_mapping_file,
    gtf_file, pfam_domains_file, interpro_entry_types_file: str
        Reference files, as for Gene.

    Returns
    -------
    list[tuple[Domain, str]]
        (domain, source transcript ID) pairs.
    """
    transcript_ids = list(superisoform.transcript_segments)
    items = []
    wanted_types = ({"DNA-binding"} if "dbi" in domain_filter else set()) | \
                   ({"PPI"} if "ppi_domain" in domain_filter else set())
    if wanted_types:
        protein_ids = {isoform["id"]: isoform["protein_id"]
                       for isoform in gene._get_gtf_index(gtf_file).get(gene.ensg_id, [])}
        for transcript_id in transcript_ids:
            for domain in Transcript.interpro_domain_objects(protein_ids.get(transcript_id), interpro_domains_file,
                                                             pfam_domains_file, interpro_entry_types_file):
                if wanted_types & set(domain.types):
                    items.append((domain, transcript_id))
    if "ppi_bs" in domain_filter:
        if ppi_binding_site_file is None:
            raise ValueError("Need ppi_binding_site_file if using PPI binding site")
        ppi_binding_sites = Transcript._get_ppi_binding_sites(ppi_binding_site_file)
        sources = choose_source_transcripts(gene.ensg_id, transcript_ids,
                                            load_uniprot_xrefs(uniprot_mapping_file),
                                            accessions=ppi_binding_sites)
        for accession, (transcript_id, _) in sorted(sources.items()):
            for domain_id, positions, source in ppi_binding_sites.get(accession, []):
                domain = Domain.from_positions_string(domain_id, positions, source)
                domain.types = ["PPI"]
                items.append((domain, transcript_id))
    return items


def map_domain(gene_id, domain, source_transcript_id, superisoform):
    """
    Map a domain's residues (in its source transcript's protein
    coordinates) to codon bases in the superisoform.

    Parameters
    ----------
    gene_id: str
    domain: Domain
    source_transcript_id: str
    superisoform: Superisoform

    Returns
    -------
    list[tuple]
        One codon (see Superisoform.residue_codons) per residue that
        falls inside the source protein; residues past its end are
        dropped with a warning.
    """
    codons = superisoform.residue_codons(source_transcript_id)
    residues = sorted({residue for part in domain.pos.parts for residue in range(part.start, part.end)})
    in_range = [residue for residue in residues if 0 <= residue < len(codons)]
    if len(in_range) < len(residues):
        warnings.warn("%s: %d of %d residues of %s lie outside %s's %d-residue protein; dropped"
                      % (gene_id, len(residues) - len(in_range), len(residues), domain.domain_id,
                         source_transcript_id, len(codons)))
    return [codons[residue] for residue in in_range]


def _merge_overlapping(mapped):
    """
    Within each classification, collapse domains sharing any codon base
    down to the largest in each overlapping group - the same greedy
    procedure as Gene.check_domain_redundancy, but comparing exact
    codon bases instead of protein start/end positions (which aren't
    comparable across transcripts). A domain of both types kept in the
    DNA-binding pass also represents its residues in the PPI pass, so
    PPI domains overlapping it are dropped rather than kept as a second
    domain on the same residues.
    """
    keeping = []
    for classification in ["DNA-binding", "PPI"]:
        kept_bases = [{base for codon in entry[2] for base in codon}
                      for entry in keeping if classification in entry[0].types]
        queue = [entry for entry in mapped if classification in entry[0].types
                 and not any({base for codon in entry[2] for base in codon} & bases for bases in kept_bases)]
        while queue:
            current = queue.pop()
            current_bases = {base for codon in current[2] for base in codon}
            overlapping = [i for i, entry in enumerate(queue)
                           if current_bases & {base for codon in entry[2] for base in codon}]
            for i in overlapping:
                if len(queue[i][2]) > len(current[2]):
                    current = queue[i]
            for i in reversed(overlapping):
                del queue[i]
            if current not in keeping:
                keeping.append(current)
    return keeping


def segment_coverage_rows(gene, superisoform, domain_filter, merge_overlapping=True,
                          ppi_binding_site_file=None, interpro_domains_file=None,
                          uniprot_mapping_file=None, gtf_file=None,
                          pfam_domains_file=None, interpro_entry_types_file=None):
    """
    Score every transcript of a gene against every one of its domains:
    coverage is the fraction of the domain's residues whose exact codon
    the transcript encodes (Superisoform.contains_codon).

    Parameters
    ----------
    gene: Gene
    superisoform: Superisoform
    domain_filter: list
        As for collect_domains.
    merge_overlapping: bool, optional
        If True (default), collapse overlapping domains within each
        classification down to the largest; if False, keep every
        domain. Exact duplicates (same domain ID and codons, e.g. one
        InterPro domain on several transcripts sharing its exons) are
        always collapsed.
    ppi_binding_site_file, interpro_domains_file, uniprot_mapping_file,
    gtf_file, pfam_domains_file, interpro_entry_types_file: str
        Reference files, as for Gene.

    Returns
    -------
    list[dict]
        One row per (transcript, domain), with keys COVERAGE_COLUMNS.
    """
    mapped = []
    seen = set()
    for domain, source in collect_domains(gene, superisoform, domain_filter, ppi_binding_site_file,
                                          interpro_domains_file, uniprot_mapping_file, gtf_file,
                                          pfam_domains_file, interpro_entry_types_file):
        codons = map_domain(gene.ensg_id, domain, source, superisoform)
        key = (domain.domain_id, frozenset(codons))
        if key in seen:
            continue
        seen.add(key)
        mapped.append((domain, source, codons))
    if merge_overlapping:
        mapped = _merge_overlapping(mapped)

    rows = []
    for transcript_id in superisoform.transcript_segments:
        for domain, source, codons in mapped:
            n_covered = sum(superisoform.contains_codon(transcript_id, codon) for codon in codons)
            rows.append({
                "gene_id": gene.ensg_id,
                "transcript_id": transcript_id,
                "domain_id": domain.domain_id,
                "domain_type": ";".join(domain.types),
                "source_transcript_id": source,
                "positions": domain.positions_str,
                "n_residues": len(codons),
                "n_covered": n_covered,
                "coverage": n_covered / len(codons) if codons else float("nan"),
            })
    return rows