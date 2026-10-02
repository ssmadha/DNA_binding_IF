import os

import pandas as pd
from Bio import SeqFeature

# Resolved relative to this file (bin/TF_ASIF/domain.py -> repo root), so the
# default works regardless of the caller's working directory (repo root for
# tests, a work/xx/<hash>/ task directory under Nextflow).
DEFAULT_DNA_BINDING_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..",
    "reference_data", "interpro_superfamily_domains_DBD.tsv")
# Pfam family pairs with a known interacting structure in 3did.
DEFAULT_PPI_PFAM_PAIRS_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..",
    "reference_data", "3did_pfam_pairs.txt")
# InterPro entry types that describe a structural unit and so may be
# called PPI. Excludes Family (whole-protein classifications, e.g. "p53
# tumour suppressor family", which would otherwise inherit PPI from any
# interacting domain inside them) and the short site types (conserved,
# active and binding sites, PTMs).
PPI_ENTRY_TYPES = {"Domain", "Homologous_superfamily", "Repeat"}
# Minimum share of both the domain and the 3did Pfam hit that their
# overlap must cover for the domain to be called PPI.
PPI_MIN_OVERLAP = 0.5

# DNA-binding classification tables, keyed by file path; read once per process.
_dna_binding_tables = {}
# Pfam families in 3did, keyed by file path; read once per process.
_ppi_pfam_families = {}


def load_ppi_pfam_families(ppi_pfam_pairs_file=DEFAULT_PPI_PFAM_PAIRS_FILE):
    """
    Load and cache the Pfam families that take part in at least one
    3did domain-domain interaction.

    Parameters
    ----------
    ppi_pfam_pairs_file: str
        Tab-separated file of interacting Pfam accession pairs, one
        pair per line (3did_pfam_pairs.txt).

    Returns
    -------
    set[str]
        Pfam accessions, e.g. "PF00046".
    """
    if ppi_pfam_pairs_file not in _ppi_pfam_families:
        families = set()
        with open(ppi_pfam_pairs_file) as handle:
            for line in handle:
                families.update(field.strip() for field in line.split("\t")[:2] if field.strip())
        _ppi_pfam_families[ppi_pfam_pairs_file] = families
    return _ppi_pfam_families[ppi_pfam_pairs_file]


def _merge_spans(spans):
    """Union of 0-based half-open (start, end) spans, as sorted disjoint spans."""
    merged = []
    for start, end in sorted(spans):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def _shared_length(spans_a, spans_b):
    """Residues shared by two lists of spans, each disjoint within itself."""
    return sum(max(0, min(end_a, end_b) - max(start_a, start_b))
               for start_a, end_a in spans_a for start_b, end_b in spans_b)


class Domain:
    """
    Domain object
    """

    def __init__(self, domain_id, positions, source, ppi_regions=None, entry_type=None):
        """
        Constructor

        Parameters
        ----------
        domain_id: domain identifier - an InterPro accession for
            InterPro-sourced domains, or the synthetic domain_id of a
            PPI binding site.
        positions: Bio.SeqFeature.Location
            This domain's position(s) - a FeatureLocation for a single
            contiguous span, or a CompoundLocation for a discontinuous
            one (e.g. a PPI contact interface made of scattered
            residues). start/end are derived from this.
        source: source
        ppi_regions: dict[str, list[tuple[int, int]]], optional
            0-based half-open spans of 3did Pfam hits on the same
            protein, keyed by Pfam accession, used by
            determine_protein_interaction. None (the default) means no
            PPI evidence, so the domain is never PPI.
        entry_type: str, optional
            InterPro entry type of domain_id (e.g. "Domain", "Family"),
            used by determine_protein_interaction.
        """
        self.domain_id = domain_id
        self.pos = positions
        self.start = positions.start
        self.end = positions.end
        self.source = source
        self.ppi_regions = ppi_regions
        self.entry_type = entry_type
        self.types = self.determine_types()

    @classmethod
    def from_positions_string(cls, domain_id, positions_str, source, ppi_regions=None, entry_type=None):
        """
        Build a Domain from a run-length-encoded positions string, as
        stored in the shared reference_data/ domain TSV schema
        (protein_id, domain_id, positions, source) used by both
        Homo_sapiens.GRCh38.interpro_domains.tsv.gz and
        ppi_binding_sites.tsv.

        Parameters
        ----------
        domain_id: domain identifier - an InterPro accession, or for a
            PPI binding site, its synthetic domain_id.
        positions_str: str
            Comma-separated list of "start-end" segments or lone
            positions, e.g. "104-352" (one contiguous span) or
            "23-37,39-43,61,94,96,99-110" (a discontinuous site).
            Positions are 1-based and inclusive (UniProt/InterPro
            residue numbering); they are stored as 0-based half-open
            locations, so "104-352" becomes FeatureLocation(103, 352)
            and "61" becomes FeatureLocation(60, 61).
        source: source
        ppi_regions, entry_type: optional
            As for __init__.

        Returns
        -------
        Domain
        """
        locations = []
        for segment in positions_str.split(","):
            if "-" in segment:
                start, end = segment.split("-")
                locations.append(SeqFeature.FeatureLocation(int(start) - 1, int(end)))
            else:
                pos = int(segment)
                locations.append(SeqFeature.FeatureLocation(pos - 1, pos))
        positions = locations[0] if len(locations) == 1 else SeqFeature.CompoundLocation(locations)
        domain = cls(domain_id, positions, source, ppi_regions=ppi_regions, entry_type=entry_type)
        domain.positions_str = positions_str
        return domain

    def determine_types(self):
        """
        Determine the types of this domain

        Returns
        -------
        list[str]
            Subset of ["DNA-binding", "PPI"] that this domain belongs
            to, based on determine_dna_binding and
            determine_protein_interaction.
        """
        types = []
        if self.determine_dna_binding():
            types.append("DNA-binding")
        if self.determine_protein_interaction():
            types.append("PPI")
        return types

    def determine_dna_binding(self, dna_binding_file=DEFAULT_DNA_BINDING_FILE):
        """
        Determine if this domain is a DNA-binding domain

        Parameters
        ----------
        dna_binding_file: str
            file containing DNA-binding domains

        Returns
        -------
        bool
            True if this domain's domain_id is listed as DNA-binding
            in dna_binding_file. False if the domain has no domain_id,
            is not sourced from SuperFamily, or is not found in the
            file.
        """
        if dna_binding_file not in _dna_binding_tables:
            _dna_binding_tables[dna_binding_file] = pd.read_csv(dna_binding_file, sep='\t', index_col=0)
        interpro_superfamily_domains_DBD = _dna_binding_tables[dna_binding_file]
        if (self.domain_id is None or self.source!="SuperFamily" or
                self.domain_id not in interpro_superfamily_domains_DBD.index):
            return False
        return interpro_superfamily_domains_DBD.loc[self.domain_id,"DNA-binding"]

    def determine_protein_interaction(self):
        """
        Determine if this domain is a protein-protein interaction
        domain: it overlaps the Pfam hits of one family with a known
        interacting structure in 3did (self.ppi_regions) by at least
        PPI_MIN_OVERLAP of both the domain's residues and those hits'
        residues. Hits of a family are pooled - only the ones touching
        this domain - so a domain spanning several repeats (e.g. a C2H2
        zinc finger superfamily entry over two single-finger hits)
        counts them together. Adapted from getTFIsoforms.ipynb, which
        needed only one shared residue. Pfam hits are only evidence
        here, never domains themselves, so a domain that is both
        DNA-binding and PPI stays one domain with both types.

        Returns
        -------
        bool
            True if this domain's InterPro entry type is in
            PPI_ENTRY_TYPES and, for some Pfam family in
            self.ppi_regions, the union of its hits touching this domain
            shares at least PPI_MIN_OVERLAP of that union's length and
            of this domain's length (summed over its parts) with this
            domain. False otherwise, including when no ppi_regions were
            given.
        """
        if not self.ppi_regions or self.entry_type not in PPI_ENTRY_TYPES:
            return False
        parts = [(int(part.start), int(part.end)) for part in self.pos.parts]
        domain_length = sum(end - start for start, end in parts)
        for hits in self.ppi_regions.values():
            touching = _merge_spans([hit for hit in hits if _shared_length([hit], parts)])
            if not touching:
                continue
            shared = _shared_length(touching, parts)
            if shared >= PPI_MIN_OVERLAP * domain_length and \
                    shared >= PPI_MIN_OVERLAP * sum(end - start for start, end in touching):
                return True
        return False

    def __repr__(self):
        """
        Build a human-readable representation of this domain.

        Returns
        -------
        str
            String with this domain's domain_id, position, and
            types.
        """
        return "Domain ID %s at %s of types %s" % (self.domain_id, self.pos, self.types)

