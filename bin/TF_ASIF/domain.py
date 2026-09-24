import pandas as pd
from Bio import SeqFeature


class Domain:
    """
    Domain object
    """

    def __init__(self, interpro_id, start, end, source, pos=None, **kwargs):
        """
        Constructor

        Parameters
        ----------
        interpro_id: interpro id
        start: start position
        end: end position
        source: source
        pos: position
        """
        self.interpro_id = interpro_id
        self.start = start
        self.end = end
        self.source = source
        self.types = self.determine_types()
        if pos is None:
            self.pos = SeqFeature.SeqFeature(SeqFeature.FeatureLocation(start, end))
        else:
            self.pos = pos

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

    def determine_dna_binding(self, dna_binding_file="/mnt/data/storage/WPI/Korkin_Lab/DNA_Binding_IF/reference_data/interpro_superfamily_domains_DBD.tsv"):
        """
        Determine if this domain is a DNA-binding domain

        Parameters
        ----------
        dna_binding_file: str
            file containing DNA-binding domains

        Returns
        -------
        bool
            True if this domain's interpro_id is listed as DNA-binding
            in dna_binding_file. False if the domain has no interpro_id,
            is not sourced from SuperFamily, or is not found in the
            file.
        """
        interpro_superfamily_domains_DBD = pd.read_csv(dna_binding_file, sep='\t', index_col=0)
        if (self.interpro_id is None or self.source!="SuperFamily" or
                self.interpro_id not in interpro_superfamily_domains_DBD.index):
            return False
        return interpro_superfamily_domains_DBD.loc[self.interpro_id,"DNA-binding"]

    def determine_protein_interaction(self):
        """
        Determine if this domain is a protein-protein interaction
        domain. Not yet implemented.

        Returns
        -------
        bool
            Always False.
        """
        return False

    def __repr__(self):
        """
        Build a human-readable representation of this domain.

        Returns
        -------
        str
            String with this domain's interpro_id, position, and
            types.
        """
        return "Interpro ID %s at %s of types %s" % (self.interpro_id, self.pos, self.types)

