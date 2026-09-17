import ensembl_rest
import pandas as pd
from Bio import SeqFeature, Align

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

    def __init__(self, gene, enst_id: str, ensp_id: str, domain_types: list, binding_site_df: pd.DataFrame,
                 idmapping_df: pd.DataFrame):
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
        idmapping_df: pd.DataFrame
            UniProt ID mapping data frame, used to resolve this
            transcript's UniProt and RefSeq IDs.
        """
        self.gene = gene
        self.enst_id = enst_id
        self.ensp_id = ensp_id
        # print("Ensembl ID: " + self.enst_id)
        self.uniprot_id = self.get_uniprot_id(idmapping_df=idmapping_df)
        # print("UniProt ID: " + str(self.uniprot_id))
        self.refseq_id = self.get_refseq_id(idmapping_df=idmapping_df)
        if self.uniprot_id is None or self.refseq_id is None:
            return
        # print("RefSeq ID: " + str(self.refseq_id))
        #self.prot_seq = self.download_sequence()
        #print("Sequence: " + self.prot_seq)
        self.domains = self.download_domains(domain_types=domain_types)
        #print("Domains: " + str(self.domains))
        #print(len(self.domains))

    def download_sequence(self, ensp_id=None):
        """
        Download the protein sequence of this transcript

        Parameters
        ----------
        ensp_id: Ensembl Protein ID

        Returns
        -------
        str
            Protein sequence for ensp_id, from the Ensembl REST API.
        """
        if ensp_id is None:
            ensp_id = self.ensp_id
        seq = ensembl_rest.sequence_id(ensp_id)["seq"]
        return seq

    def get_uniprot_id(self, idmapping_df: pd.DataFrame, ensp_id=None):
        """
        Identify the uniprot id of this transcript

        Parameters
        ----------
        idmapping_df: pd.DataFrame
            UniProt ID mapping data frame to search for a row mapping
            ensp_id to a UniProt ID.
        ensp_id: Ensembl Protein ID

        Returns
        -------
        str or None
            The first matching UniProt ID, or None if ensp_id is not
            found in idmapping_df.
        """
        if ensp_id is None:
            ensp_id = self.ensp_id

        uniprot_ids = idmapping_df.loc[(idmapping_df[2].str.contains(ensp_id)) &
                                       (idmapping_df[1]=="Ensembl_PRO"), 0].tolist()
        if len(uniprot_ids)>0:
            return uniprot_ids[0]
        else:
            return None

    def get_refseq_id(self, idmapping_df: pd.DataFrame, uniprot_id: str | None = None):
        """
        Identify the refseq id of this transcript

        Parameters
        ----------
        idmapping_df: pd.DataFrame
            UniProt ID mapping data frame to search for a row mapping
            uniprot_id to a RefSeq protein ID.
        uniprot_id: UniProt ID

        Returns
        -------
        str or None
            The first matching RefSeq protein ID (starting with "NP_"),
            or None if uniprot_id is not found in idmapping_df or is
            unavailable.
        """
        if uniprot_id is None:
            if self.uniprot_id is not None:
                uniprot_id: str = self.uniprot_id
            else:
                return None
        refseq_ids = idmapping_df.loc[(idmapping_df[0].str.contains(uniprot_id)) &
                                      (idmapping_df[1]=="RefSeq") &
                                      (idmapping_df[2].str.startswith("NP_")), 2].tolist()
        if len(refseq_ids)>0:
            return refseq_ids[0]
        else:
            return None

    def download_domains(self, ensp_id=None, domain_types=None, binding_site_df: pd.DataFrame | None = None):
        """
        Download this transcript's domains from Ensembl and/or the
        PPI binding site data, depending on domain_types.

        Parameters
        ----------
        ensp_id : str, optional
            Ensembl Protein ID. Defaults to self.ensp_id.
        domain_types : list, optional
            Domain types to include. "ppi_domain" or "dbi" fetch
            SuperFamily domains from the Ensembl REST API; "ppi_bs"
            fetches PPI binding sites via yue_ppi_locations. Defaults
            to ["ppi_domain", "dbi"].
        binding_site_df : pd.DataFrame, optional
            Data frame of PPI binding sites, required when "ppi_bs"
            is in domain_types.

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
        if "ppi_domain" in domain_types or "dbi" in domain_types:
            results = ensembl_rest.overlap_translation(ensp_id,
                                                   type="domain")
            domains += [Domain(interpro_id=res["interpro"], source = res["type"], **res)
                        for res in results if res["type"] == "SuperFamily"]
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
            parent gene, UniProt ID, and (if resolved) RefSeq ID and
            exon locations.
        """
        return_string = "Transcript Ensembl ID: {}".format(self.enst_id)
        return_string += "\n Part of Gene: {}".format(self.gene.ensg_id)
        return_string += "\n Uniprot ID: {}".format(self.uniprot_id)
        if self.refseq_id is not None:
            return_string += "\n RefSeq ID: {}".format(self.refseq_id)
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
            parent gene, UniProt ID, and (if resolved) RefSeq ID and
            exon locations.
        """
        return_string = "Transcript Ensembl ID: {}".format(self.enst_id)
        return_string += "\n Part of Gene: {}".format(self.gene.ensg_id)
        return_string += "\n Uniprot ID: {}".format(self.uniprot_id)
        if self.refseq_id is not None:
            return_string += "\n RefSeq ID: {}".format(self.refseq_id)
            return_string += "\n RNA Exons at: {}".format(self.exons_rna)
            return_string += "\n Protein Exons at: {}".format(self.exons_prot)
        return return_string
