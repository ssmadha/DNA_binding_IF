import pandas as pd
from Bio import SeqFeature


class Domain:
    """
    Domain object
    """

    def __init__(self, interpro_id, positions, source):
        """
        Constructor

        Parameters
        ----------
        interpro_id: interpro id
        positions: Bio.SeqFeature.Location
            This domain's position(s) - a FeatureLocation for a single
            contiguous span, or a CompoundLocation for a discontinuous
            one (e.g. a PPI contact interface made of scattered
            residues). start/end are derived from this.
        source: source
        """
        self.interpro_id = interpro_id
        self.pos = positions
        self.start = positions.start
        self.end = positions.end
        self.source = source
        self.types = self.determine_types()

    @classmethod
    def from_positions_string(cls, interpro_id, positions_str, source):
        """
        Build a Domain from a run-length-encoded positions string, as
        stored in the shared reference_data/ domain TSV schema
        (protein_id, domain_id, positions, source) used by both
        Homo_sapiens.GRCh38.interpro_domains.tsv.gz and
        ppi_binding_sites.tsv.

        Parameters
        ----------
        interpro_id: interpro id (or, for a PPI binding site, its
            synthetic domain_id).
        positions_str: str
            Comma-separated list of "start-end" segments or lone
            positions, e.g. "104-352" (one contiguous span) or
            "23-37,39-43,61,94,96,99-110" (a discontinuous site).
        source: source

        Returns
        -------
        Domain
        """
        locations = []
        for segment in positions_str.split(","):
            if "-" in segment:
                start, end = segment.split("-")
                locations.append(SeqFeature.FeatureLocation(int(start), int(end)))
            else:
                pos = int(segment)
                locations.append(SeqFeature.FeatureLocation(pos, pos))
        positions = locations[0] if len(locations) == 1 else SeqFeature.CompoundLocation(locations)
        return cls(interpro_id, positions, source)

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

