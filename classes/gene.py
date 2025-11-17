
import ensembl_rest
import mygene


class Domain:
    """
    Domain object
    """

    def __init__(self, interpro_id, start, end, source, classification = None):
        self.interpro_id = interpro_id
        self.classification = classification
        self.start = start
        self.end = end
        self.source = source

class Transcript:
    """
    Transcript object
    """
    seq = None
    uniprot_id = None
    domains = None

    def __init__(self, ensembl_id: str):
        self.ensembl_id = ensembl_id
        self.seq = self.download_sequence()

    def download_sequence(self, ensembl_id=None):
        if ensembl_id==None:
            ensembl_id = self.ensembl_id
        seq = ensembl_rest.sequence_id(ensembl_id)
        return seq

    def download_domains(self):
        results = ensembl_rest.overlap_translation(self.ensembl_id, 
                                                   type="domain")
        self.domains = [res for res in results]

class Gene:
    """
    Gene object
    """
    seq = None
    uniprot_id = None
    transcripts = None
    
    def __init__(self, ensembl_id: str, biotype_filter=None):
        """
        Initialize a Gene object based on Ensembl ID

        params:

        """
        if biotype_filter is None:
            biotype_filter = ['protein_coding']
        self.ensembl_id = ensembl_id
        self.transcripts = self.download_transcripts(self.ensembl_id, biotype_filter)
        self.uniprot_id = self.download_uniprot_id()

    def download_transcripts(self, ensembl_gene_id=None, biotype_filter=None):
        if biotype_filter is None:
            biotype_filter = ['protein_coding']
        if ensembl_gene_id is None:
            ensembl_gene_id = self.ensembl_id
        try:
            isoforms = ensembl_rest.lookup(ensembl_gene_id,
                                           params={'multiple_sequences': True,
                                                   'type': 'protein',
                                                   'expand': True
                                                   }
                                           )
        except ensembl_rest.HTTPError as err:
            error_code = err.response.status_code
            error_message= err.response.json()['error']
            if(error_code==400) and ("not found" in error_message):
                print(ensembl_gene_id + " not found")
            elif(error_code==400) and ("No sequences returned" in error_message):
                print(ensembl_gene_id + " no protein sequences found")
            else:
                raise
            return None
        transcripts = []
        for isoform in isoforms["Transcript"]:
            if isoform['biotype'] not in biotype_filter:
                continue
            transcript = Transcript(isoform['id'])
            transcripts.append(transcript)
        return transcripts

    def download_uniprot_id(self, ensembl_id = None):
        if ensembl_id is None:
            ensembl_id = self.ensembl_id
        mg = mygene.MyGeneInfo()
        get_gene_result = mg.getgene(ensembl_id)
        print(get_gene_result)
        uniprot_id = None
        if get_gene_result is not None and "uniprot" in get_gene_result:
            uniprot_id = get_gene_result["uniprot"]
        print("No UniProt ID found for " + ensembl_id)
        return uniprot_id