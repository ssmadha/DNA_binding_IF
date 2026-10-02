import gzip
import warnings

import pandas as pd
from Bio import Align, SeqIO

from .domain import Domain, load_ppi_pfam_families


class Transcript:
    """
    Transcript object
    """
    uniprot_id = None
    domains = []
    rna_seq = None
    exons_rna = None
    prot_seq = None
    exons_prot = None
    _cds_sequences = None
    _uniprot_mapping = None
    _interpro_domains = None
    _ppi_binding_sites = None
    _ppi_regions = None
    _interpro_entry_types = None

    def __init__(self, gene, enst_id: str, ensp_id: str, domain_types: list,
                 cds_fasta_file: str = None, uniprot_mapping_file: str = None,
                 interpro_domains_file: str = None, ppi_binding_site_file: str = None,
                 pfam_domains_file: str = None, interpro_entry_types_file: str = None):
        """
        Constructor

        Parameters
        ----------
        gene: parent gene object
        enst_id: Ensembl Transcript ID
        ensp_id: Ensembl Protein ID
        domain_types: list of domain types to use
        cds_fasta_file: str, optional
            Path to a (optionally gzipped) Ensembl "cds.all.fa" FASTA
            file, used by download_sequence for local protein sequence
            lookup.
        uniprot_mapping_file: str, optional
            Path to a (optionally gzipped) Ensembl protein-to-UniProt
            xref TSV file (e.g. "*.uniprot.tsv.gz"), used by
            get_uniprot_id to resolve this transcript's UniProt ID.
        interpro_domains_file: str, optional
            Path to a (optionally gzipped) local InterPro protein
            domain TSV file (protein_id, domain_id, positions,
            source), used by download_domains to resolve this
            transcript's SuperFamily/InterPro domains.
        ppi_binding_site_file: str, optional
            Path to a tab-separated PPI binding site file (same
            protein_id/domain_id/positions/source schema as
            interpro_domains_file, keyed by UniProt ID instead of
            Ensembl Protein ID), used by yue_ppi_locations when
            "ppi_bs" is in domain_types.
        pfam_domains_file, interpro_entry_types_file: str, optional
            Evidence for classifying InterPro domains as PPI; see
            interpro_domain_objects.
        """
        self.gene = gene
        self.enst_id = enst_id
        self.ensp_id = ensp_id
        self.cds_fasta_file = cds_fasta_file
        self.uniprot_mapping_file = uniprot_mapping_file
        self.interpro_domains_file = interpro_domains_file
        self.ppi_binding_site_file = ppi_binding_site_file
        self.pfam_domains_file = pfam_domains_file
        self.interpro_entry_types_file = interpro_entry_types_file
        # print("Ensembl ID: " + self.enst_id)
        self.uniprot_id = self.get_uniprot_id()
        # print("UniProt ID: " + str(self.uniprot_id))
        if self.uniprot_id is None:
            return
        #self.prot_seq = self.download_sequence()
        #print("Sequence: " + self.prot_seq)
        self.domains = self.download_domains(domain_types=domain_types)
        #print("Domains: " + str(self.domains))
        #print(len(self.domains))

    def download_sequence(self, enst_id=None, cds_fasta_file=None):
        """
        Look up the protein sequence of this transcript by translating
        its CDS sequence from a local Ensembl CDS FASTA file, instead
        of calling the Ensembl REST API.

        Parameters
        ----------
        enst_id: Ensembl Transcript ID
        cds_fasta_file: str, optional
            Path to a (optionally gzipped) Ensembl "cds.all.fa" FASTA
            file to look up enst_id's CDS sequence in. Defaults to
            self.cds_fasta_file.

        Returns
        -------
        str or None
            Protein sequence translated from enst_id's CDS sequence,
            or None if the CDS is partial (its length is not a
            multiple of 3) and was skipped.
        """
        if enst_id is None:
            enst_id = self.enst_id
        if cds_fasta_file is None:
            cds_fasta_file = self.cds_fasta_file
        cds_sequences = self._get_cds_sequences(cds_fasta_file)
        cds_seq = cds_sequences[enst_id]
        if len(cds_seq) % 3 != 0:
            warnings.warn("Skipping " + enst_id + ": partial CDS (length "
                          + str(len(cds_seq)) + " is not a multiple of 3)")
            return None
        return str(cds_seq.translate())

    @classmethod
    def _get_cds_sequences(cls, cds_fasta_file):
        """
        Load and cache CDS sequences from an Ensembl CDS FASTA file,
        keyed by version-less Ensembl Transcript ID. The file is only
        parsed once per process; subsequent calls reuse the cache.

        Parameters
        ----------
        cds_fasta_file: str
            Path to a (optionally gzipped) Ensembl "cds.all.fa" FASTA
            file.

        Returns
        -------
        dict[str, Bio.Seq.Seq]
            CDS sequences keyed by version-less Ensembl Transcript ID.
        """
        if cls._cds_sequences is None:
            opener = gzip.open if cds_fasta_file.endswith(".gz") else open
            with opener(cds_fasta_file, "rt") as handle:
                cls._cds_sequences = {record.id.split(".")[0]: record.seq
                                      for record in SeqIO.parse(handle, "fasta")}
        return cls._cds_sequences

    def get_uniprot_id(self, ensp_id=None, uniprot_mapping_file=None):
        """
        Identify the uniprot id of this transcript by looking up its
        Ensembl Protein ID in a local Ensembl protein-to-UniProt xref
        TSV file.

        Parameters
        ----------
        ensp_id: Ensembl Protein ID
        uniprot_mapping_file: str, optional
            Path to a (optionally gzipped) Ensembl protein-to-UniProt
            xref TSV file (e.g. "*.uniprot.tsv.gz"). Defaults to
            self.uniprot_mapping_file.

        Returns
        -------
        str or None
            The UniProt ID mapped to ensp_id, preferring the isoform-
            specific accession (e.g. "P41235-5") when one is known,
            then a reviewed SWISSPROT entry, then an unreviewed
            SPTREMBL one. None if ensp_id has no UniProt mapping.
        """
        if ensp_id is None:
            ensp_id = self.ensp_id
        if uniprot_mapping_file is None:
            uniprot_mapping_file = self.uniprot_mapping_file
        uniprot_mapping = self._get_uniprot_mapping(uniprot_mapping_file)
        return uniprot_mapping.get(ensp_id)

    @classmethod
    def _get_uniprot_mapping(cls, uniprot_mapping_file):
        """
        Load and cache a mapping of Ensembl Protein ID to UniProt ID
        from an Ensembl protein-to-UniProt xref TSV file, preferring
        the isoform-specific ("Uniprot_isoform") accession for a given
        protein over the reviewed SWISSPROT accession, which in turn
        is preferred over the unreviewed SPTREMBL one. The file is
        only parsed once per process; subsequent calls reuse the
        cache.

        Parameters
        ----------
        uniprot_mapping_file: str
            Path to a (optionally gzipped) Ensembl protein-to-UniProt
            xref TSV file (e.g. "*.uniprot.tsv.gz").

        Returns
        -------
        dict[str, str]
            UniProt IDs keyed by Ensembl Protein ID.
        """
        if cls._uniprot_mapping is None:
            mapping_df = pd.read_csv(uniprot_mapping_file, sep='\t')
            db_priority = {"Uniprot_isoform": 0, "Uniprot/SWISSPROT": 1, "Uniprot/SPTREMBL": 2}
            mapping_df = mapping_df[mapping_df["db_name"].isin(db_priority)]
            mapping_df = mapping_df.sort_values(
                by="db_name", key=lambda col: col.map(db_priority))
            mapping_df = mapping_df.drop_duplicates(subset="protein_stable_id", keep="first")
            cls._uniprot_mapping = dict(zip(mapping_df["protein_stable_id"], mapping_df["xref"]))
        return cls._uniprot_mapping

    @classmethod
    def _load_domain_reference_file(cls, file_path):
        """
        Load a shared-schema domain reference TSV file (columns:
        protein_id, domain_id, positions, source - see
        Domain.from_positions_string for the "positions" format),
        keyed by protein_id. Used for both
        Homo_sapiens.GRCh38.interpro_domains.tsv.gz (keyed by Ensembl
        Protein ID) and ppi_binding_sites.tsv (keyed by UniProt ID).

        Parameters
        ----------
        file_path: str
            Path to a (optionally gzipped) domain reference TSV file.

        Returns
        -------
        dict[str, list[tuple[str, str, str]]]
            For each protein_id, a list of (domain_id, positions,
            source) tuples.
        """
        domains_df = pd.read_csv(file_path, sep='\t')
        domains = {}
        for protein_id, domain_id, positions, source in zip(
                domains_df["protein_id"], domains_df["domain_id"],
                domains_df["positions"], domains_df["source"]):
            domains.setdefault(protein_id, []).append((domain_id, positions, source))
        return domains

    @classmethod
    def _get_interpro_domains(cls, interpro_domains_file):
        """
        Load and cache a mapping of Ensembl Protein ID to its InterPro
        protein domains from a local InterPro protein domain TSV file
        (an Ensembl BioMart "Protein Domains and Families" InterPro
        export). The file is only parsed once per process; subsequent
        calls reuse the cache.

        Parameters
        ----------
        interpro_domains_file: str
            Path to a (optionally gzipped) local InterPro protein
            domain TSV file.

        Returns
        -------
        dict[str, list[tuple[str, str, str]]]
            See _load_domain_reference_file.
        """
        if cls._interpro_domains is None:
            cls._interpro_domains = cls._load_domain_reference_file(interpro_domains_file)
        return cls._interpro_domains

    @classmethod
    def _get_ppi_binding_sites(cls, ppi_binding_site_file):
        """
        Load and cache a mapping of UniProt ID to PPI binding site
        domains from a local PPI binding site TSV file. The file is
        only parsed once per process; subsequent calls reuse the
        cache.

        Parameters
        ----------
        ppi_binding_site_file: str
            Path to a tab-separated PPI binding site file.

        Returns
        -------
        dict[str, list[tuple[str, str, str]]]
            See _load_domain_reference_file.
        """
        if cls._ppi_binding_sites is None:
            cls._ppi_binding_sites = cls._load_domain_reference_file(ppi_binding_site_file)
        return cls._ppi_binding_sites

    @classmethod
    def _get_ppi_regions(cls, pfam_domains_file):
        """
        Load and cache, for each Ensembl Protein ID, the spans of its
        Pfam hits whose family is in 3did (see
        domain.load_ppi_pfam_families), grouped by Pfam family. The file
        is only parsed once per process; subsequent calls reuse the
        cache.

        Parameters
        ----------
        pfam_domains_file: str
            Path to a (optionally gzipped) Pfam hit TSV file in the
            shared domain schema (Homo_sapiens.GRCh38.pfam_domains.tsv.gz).

        Returns
        -------
        dict[str, dict[str, list[tuple[int, int]]]]
            Per protein, 0-based half-open (start, end) hit spans per
            Pfam accession.
        """
        if cls._ppi_regions is None:
            families = load_ppi_pfam_families()
            cls._ppi_regions = {}
            for protein_id, hits in cls._load_domain_reference_file(pfam_domains_file).items():
                regions = {}
                for pfam_id, positions, _ in hits:
                    if pfam_id not in families:
                        continue
                    location = Domain.from_positions_string(pfam_id, positions, "Pfam").pos
                    regions.setdefault(pfam_id, []).extend(
                        (int(part.start), int(part.end)) for part in location.parts)
                if regions:
                    cls._ppi_regions[protein_id] = regions
        return cls._ppi_regions

    @classmethod
    def _get_interpro_entry_types(cls, interpro_entry_types_file):
        """
        Load and cache the type (Domain, Family, Homologous_superfamily,
        ...) of every InterPro entry from an InterPro entry.list file
        (columns ENTRY_AC, ENTRY_TYPE, ENTRY_NAME).

        Parameters
        ----------
        interpro_entry_types_file: str
            Path to an InterPro entry.list file.

        Returns
        -------
        dict[str, str]
            InterPro accession -> entry type.
        """
        if cls._interpro_entry_types is None:
            entries_df = pd.read_csv(interpro_entry_types_file, sep='\t')
            cls._interpro_entry_types = dict(zip(entries_df["ENTRY_AC"], entries_df["ENTRY_TYPE"]))
        return cls._interpro_entry_types

    @classmethod
    def interpro_domain_objects(cls, ensp_id, interpro_domains_file,
                                pfam_domains_file=None, interpro_entry_types_file=None):
        """
        Build the classified InterPro domains of one Ensembl protein.

        Parameters
        ----------
        ensp_id: str
            Ensembl Protein ID.
        interpro_domains_file: str
            See _get_interpro_domains.
        pfam_domains_file, interpro_entry_types_file: str, optional
            See _get_ppi_regions and _get_interpro_entry_types. Both are
            needed to classify domains as PPI; if neither is given, no
            domain is PPI.

        Returns
        -------
        list[Domain]
        """
        if (pfam_domains_file is None) != (interpro_entry_types_file is None):
            raise ValueError("PPI domain classification needs both pfam_domains_file "
                             "and interpro_entry_types_file")
        ppi_regions, entry_types = None, {}
        if pfam_domains_file is not None:
            ppi_regions = cls._get_ppi_regions(pfam_domains_file).get(ensp_id)
            entry_types = cls._get_interpro_entry_types(interpro_entry_types_file)
        return [Domain.from_positions_string(domain_id, positions, source, ppi_regions=ppi_regions,
                                             entry_type=entry_types.get(domain_id))
                for domain_id, positions, source in cls._get_interpro_domains(interpro_domains_file).get(ensp_id, [])]

    def download_domains(self, ensp_id=None, domain_types=None,
                         interpro_domains_file=None, ppi_binding_site_file=None):
        """
        Look up this transcript's domains in a local InterPro protein
        domain file and/or the PPI binding site data, depending on
        domain_types.

        Parameters
        ----------
        ensp_id : str, optional
            Ensembl Protein ID. Defaults to self.ensp_id.
        domain_types : list, optional
            Domain types to include. "dbi" and "ppi_domain" look up
            SuperFamily/InterPro domains in a local InterPro protein
            domain file, keeping those classified DNA-binding and PPI
            respectively; "ppi_bs" fetches PPI binding sites via
            yue_ppi_locations. Defaults to ["ppi_domain", "dbi"].
        interpro_domains_file : str, optional
            Path to a (optionally gzipped) local InterPro protein
            domain TSV file. Defaults to self.interpro_domains_file.
        ppi_binding_site_file : str, optional
            Path to a tab-separated PPI binding site file, required
            when "ppi_bs" is in domain_types. Defaults to
            self.ppi_binding_site_file.

        Returns
        -------
        list[Domain]
            Domains found for this transcript.
        """
        if domain_types is None:
            domain_types = ["ppi_domain", "dbi"]
        domains = []
        if ensp_id is None:
            ensp_id = self.ensp_id
        if interpro_domains_file is None:
            interpro_domains_file = self.interpro_domains_file
        if ppi_binding_site_file is None:
            ppi_binding_site_file = self.ppi_binding_site_file
        wanted_types = ({"DNA-binding"} if "dbi" in domain_types else set()) | \
                       ({"PPI"} if "ppi_domain" in domain_types else set())
        if wanted_types:
            domains += [domain for domain in self.interpro_domain_objects(
                            ensp_id, interpro_domains_file, self.pfam_domains_file, self.interpro_entry_types_file)
                        if wanted_types & set(domain.types)]
        if "ppi_bs" in domain_types:
            if ppi_binding_site_file is None:
                raise ValueError("Need ppi_binding_site_file if using PPI binding site")
            domains += self.yue_ppi_locations(ppi_binding_site_file=ppi_binding_site_file)
        return domains

    def yue_ppi_locations(self, ppi_binding_site_file=None) -> list[Domain]:
        """
        Identify PPI locations

        Parameters
        ----------
        ppi_binding_site_file : str, optional
            Path to a tab-separated PPI binding site file. Defaults to
            self.ppi_binding_site_file.

        Returns
        -------
        list[Domain]
            List of PPI Domains
        """
        if self.uniprot_id is None:
            return []
        if ppi_binding_site_file is None:
            ppi_binding_site_file = self.ppi_binding_site_file
        uniprot_id = self.uniprot_id.split("-")[0]
        ppi_binding_sites = self._get_ppi_binding_sites(ppi_binding_site_file)
        domains = [Domain.from_positions_string(domain_id, positions, source)
                   for domain_id, positions, source in ppi_binding_sites.get(uniprot_id, [])]
        for domain in domains:
            domain.types = ["PPI"]
        return domains

    def align_to_reference(self, refmode="superisoform", alignmode="global", identical_only=False):
        """
        Globally align this transcript's protein sequence to a
        reference sequence, and return the coverage percentage of each
        superdomain over the best-covering alignment.

        Parameters
        ----------
        refmode : str, optional
            If "superisoform", aligns against self.gene.superisoform_seq.
            Otherwise, aligns against the protein sequence of the
            gene's first transcript.
        alignmode : str, optional
            Intended alignment mode for the pairwise aligner. Currently
            unused; the aligner mode is hardcoded to "global".
        identical_only : bool, optional
            If True, a reference residue only counts as covered when it
            is aligned to an identical residue; mismatches count as
            uncovered, like gaps. If False (default), any aligned
            residue counts as covered.

        Returns
        -------
        dict[str, float]
            Coverage fraction keyed by superdomain ID.
        """
        aligner = Align.PairwiseAligner()
        aligner.match_score = 10
        aligner.mismatch_score = -15
        aligner.open_insertion_score = -20
        aligner.extend_insertion_score = -20
        aligner.open_deletion_score = -25
        aligner.extend_deletion_score = 0
        aligner.mode = "global"
        if refmode=="superisoform":
            ref_seq = self.gene.superisoform_seq
        else:
            ref_seq = self.gene.transcripts[0].prot_seq
        #print(ref_seq)
        transcript_seq = self.prot_seq
        #print(transcript_seq)
        alignments = aligner.align(ref_seq, transcript_seq)
        superdomains = self.gene.superdomains
        #print([domain.pos for domain in superdomains])

        isoform_coverage_percentages = {}
        for i in range(len(alignments)):
            # One character per reference position (domain positions are in
            # reference coordinates): the transcript residue aligned there,
            # or "-" if none (or, with identical_only, if it differs).
            # Transcript residues with no reference position (insertions)
            # are skipped rather than shifting later positions.
            target_indices, query_indices = alignments[i].indices
            aligned_query = ["-"] * len(ref_seq)
            for k, j in zip(target_indices, query_indices):
                if k != -1 and j != -1 and not (identical_only and ref_seq[k] != transcript_seq[j]):
                    aligned_query[k] = transcript_seq[j]
            aligned_query = "".join(aligned_query)
            for domain in superdomains:
                domain_query = domain.pos.extract(aligned_query)
                # print(domain_query)
                # print(len(domain_query))
                # print(alignments[i].counts())
                overlap_perc = 1 - domain_query.count("-") / len(domain_query)
                if domain.domain_id not in isoform_coverage_percentages or \
                        isoform_coverage_percentages[domain.domain_id] < overlap_perc:
                    # print(alignments[i])
                    isoform_coverage_percentages[domain.domain_id] = overlap_perc
        return isoform_coverage_percentages

    def __repr__(self):
        """
        Build a human-readable representation of this transcript.

        Returns
        -------
        str
            Multi-line string with this transcript's Ensembl ID,
            parent gene, UniProt ID, and (if resolved) exon locations.
        """
        return_string = "Transcript Ensembl ID: {}".format(self.enst_id)
        return_string += "\n Part of Gene: {}".format(self.gene.ensg_id)
        return_string += "\n Uniprot ID: {}".format(self.uniprot_id)
        if self.uniprot_id is not None:
            return_string += "\n RNA Exons at: {}".format(self.exons_rna)
            return_string += "\n Protein Exons at: {}".format(self.exons_prot)
        return return_string

    def __str__(self):
        """
        Build a human-readable representation of this transcript.

        Returns
        -------
        str
            Multi-line string with this transcript's Ensembl ID,
            parent gene, UniProt ID, and (if resolved) exon locations.
        """
        return_string = "Transcript Ensembl ID: {}".format(self.enst_id)
        return_string += "\n Part of Gene: {}".format(self.gene.ensg_id)
        return_string += "\n Uniprot ID: {}".format(self.uniprot_id)
        if self.uniprot_id is not None:
            return_string += "\n RNA Exons at: {}".format(self.exons_rna)
            return_string += "\n Protein Exons at: {}".format(self.exons_prot)
        return return_string
