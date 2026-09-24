import gzip
import warnings

import pandas as pd
from Bio import SeqFeature, Align, SeqIO

from bin.TF_ASIF.domain import Domain


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

    def __init__(self, gene, enst_id: str, ensp_id: str, domain_types: list, binding_site_df: pd.DataFrame,
                 cds_fasta_file: str = None, uniprot_mapping_file: str = None, interpro_domains_file: str = None):
        """
        Constructor

        Parameters
        ----------
        gene: parent gene object
        enst_id: Ensembl Transcript ID
        ensp_id: Ensembl Protein ID
        domain_types: list of domain types to use
        binding_site_df: pd.DataFrame
            Data frame of PPI binding sites, used when "ppi_bs" is in
            domain_types.
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
            domain TSV file (protein_stable_id, interpro_id, start,
            end), used by download_domains to resolve this
            transcript's SuperFamily/InterPro domains.
        """
        self.gene = gene
        self.enst_id = enst_id
        self.ensp_id = ensp_id
        self.cds_fasta_file = cds_fasta_file
        self.uniprot_mapping_file = uniprot_mapping_file
        self.interpro_domains_file = interpro_domains_file
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
    def _get_interpro_domains(cls, interpro_domains_file):
        """
        Load and cache a mapping of Ensembl Protein ID to its InterPro
        protein domains from a local InterPro protein domain TSV file
        (columns: protein_stable_id, interpro_id, start, end; an
        Ensembl BioMart "Protein Domains and Families" InterPro
        export). The file is only parsed once per process; subsequent
        calls reuse the cache.

        Parameters
        ----------
        interpro_domains_file: str
            Path to a (optionally gzipped) local InterPro protein
            domain TSV file.

        Returns
        -------
        dict[str, list[tuple[str, int, int]]]
            For each Ensembl Protein ID, a list of (interpro_id,
            start, end) tuples.
        """
        if cls._interpro_domains is None:
            domains_df = pd.read_csv(interpro_domains_file, sep='\t')
            interpro_domains = {}
            for protein_stable_id, interpro_id, start, end in zip(
                    domains_df["protein_stable_id"], domains_df["interpro_id"],
                    domains_df["start"], domains_df["end"]):
                interpro_domains.setdefault(protein_stable_id, []).append((interpro_id, int(start), int(end)))
            cls._interpro_domains = interpro_domains
        return cls._interpro_domains

    def download_domains(self, ensp_id=None, domain_types=None, binding_site_df: pd.DataFrame | None = None,
                         interpro_domains_file=None):
        """
        Look up this transcript's domains in a local InterPro protein
        domain file and/or the PPI binding site data, depending on
        domain_types.

        Parameters
        ----------
        ensp_id : str, optional
            Ensembl Protein ID. Defaults to self.ensp_id.
        domain_types : list, optional
            Domain types to include. "ppi_domain" or "dbi" look up
            SuperFamily/InterPro domains in a local InterPro protein
            domain file; "ppi_bs" fetches PPI binding sites via
            yue_ppi_locations. Defaults to ["ppi_domain", "dbi"].
        binding_site_df : pd.DataFrame, optional
            Data frame of PPI binding sites, required when "ppi_bs"
            is in domain_types.
        interpro_domains_file : str, optional
            Path to a (optionally gzipped) local InterPro protein
            domain TSV file. Defaults to self.interpro_domains_file.

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
        if "ppi_domain" in domain_types or "dbi" in domain_types:
            interpro_domains = self._get_interpro_domains(interpro_domains_file)
            domains += [Domain(interpro_id=interpro_id, source="SuperFamily", start=start, end=end)
                        for interpro_id, start, end in interpro_domains.get(ensp_id, [])]
        if "ppi_bs" in domain_types:
            if binding_site_df is None :
                raise "Need binding site df if using PPI binding site"
            else:
                domains += self.yue_ppi_locations(binding_site_df=binding_site_df)
        return domains

    def yue_ppi_locations(self, binding_site_df: pd.DataFrame) -> list[Domain]:
        """
        Identify PPI locations

        Parameters
        ----------
        binding_site_df : pd.DataFrame
            Data frame containing binding sites
        Returns
        -------
        list[Domain]
            List of PPI Domains
        """
        if self.uniprot_id is None:
            return []
        domains = []
        uniprot_id = self.uniprot_id.split("-")[0]
        for index_number in binding_site_df.index[binding_site_df["UniProt"] == uniprot_id]:
            binding_site_id = binding_site_df.loc[index_number, "ID"]
            binding_site_source = binding_site_df.loc[index_number, "Source"]
            binding_site = binding_site_df.loc[index_number, "Binding_Site"]
            locs = [SeqFeature.FeatureLocation(int(loc.split(", ")[0]), int(loc.split(", ")[-1]))
                    for loc in binding_site[1:-1].split(", ")]
            if len(locs) == 1:
                domains.append(Domain(interpro_id=binding_site_id, source=binding_site_source, start=locs[0].start, end=locs[0].end, pos=locs[0]))
            else:
                domains.append(Domain(interpro_id=binding_site_id, source=binding_site_source, start=locs[0].start, end=locs[-1].end,
                                      pos=SeqFeature.CompoundLocation(locs)))
        for domain in domains:
            domain.types=["PPI"]
        return domains

    def align_to_reference(self, refmode="superisoform", alignmode="global"):
        """
        Globally align this transcript's protein sequence to a
        reference sequence, and print the coverage percentage of each
        superdomain that overlaps the best-covering alignment.

        Parameters
        ----------
        refmode : str, optional
            If "superisoform", aligns against self.gene.superisoform_seq.
            Otherwise, aligns against the protein sequence of the
            gene's first transcript.
        alignmode : str, optional
            Intended alignment mode for the pairwise aligner. Currently
            unused; the aligner mode is hardcoded to "global".
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
            for domain in superdomains:
                domain_query = domain.pos.extract("".join([alignments[i].query[j] if j!=-1 else "-" for j in alignments[i].indices[1]]))
                # print(domain_query)
                # print(len(domain_query))
                # print(alignments[i].counts())
                overlap_perc = 1 - domain_query.count("-") / len(domain_query)
                if domain.interpro_id not in isoform_coverage_percentages or \
                        isoform_coverage_percentages[domain.interpro_id] < overlap_perc:
                    # print(alignments[i])
                    isoform_coverage_percentages[domain.interpro_id] = overlap_perc
        print(self.gene.ensg_id)
        print(self.enst_id)
        print(len(isoform_coverage_percentages))
        print(list(isoform_coverage_percentages.values()))

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
