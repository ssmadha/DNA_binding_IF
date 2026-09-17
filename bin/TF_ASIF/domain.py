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
        """
        types = []
        if self.determine_dna_binding():
            types.append("DNA-binding")
        if self.determine_protein_interaction():
            types.append("PPI")
        return types

    def determine_dna_binding(self, dna_binding_file="/mnt/data/storage/WPI/Korkin_Lab/DNA_Binding_IF/interpro_superfamily_domains_DBD.tsv"):
        """
        Determine if this domain is a DNA-binding domain

        Parameters
        ----------
        dna_binding_file: str
            file containing DNA-binding domains
        """
        interpro_superfamily_domains_DBD = pd.read_csv(dna_binding_file, sep='\t', index_col=0)
        if (self.interpro_id is None or self.source!="SuperFamily" or
                self.interpro_id not in interpro_superfamily_domains_DBD.index):
            return False
        return interpro_superfamily_domains_DBD.loc[self.interpro_id,"DNA-binding"]

    def determine_protein_interaction(self):
        """

        Returns
        -------

        """
        return False

    def __repr__(self):
        return "Interpro ID %s at %s of types %s" % (self.interpro_id, self.pos, self.types)

