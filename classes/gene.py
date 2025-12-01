
import ensembl_rest
import mygene
import re


class Domain:
    """
    Domain object
    """

    def __init__(self, interpro_id, start, end, source, **kwargs):
        self.interpro_id = interpro_id
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

    def __init__(self, enst_id: str, ensp_id: str):
        self.enst_id = enst_id
        self.ensp_id = ensp_id
        self.seq = self.download_sequence()
        self.domains = self.download_domains()

    def download_sequence(self, enst_id=None):
        if enst_id==None:
            enst_id = self.enst_id
        seq = ensembl_rest.sequence_id(enst_id)
        return seq

    def download_domains(self, ensp_id=None):
        if ensp_id is None:
            ensp_id = self.ensp_id
        results = ensembl_rest.overlap_translation(ensp_id,
                                                   type="domain")
        #print(results)
        return [Domain(interpro_id=res["interpro"], source = res["type"], **res) for res in results]

class Gene:
    """
    Gene object
    """
    seq = None
    uniprot_id = None
    refseq_id = None
    transcripts = None
    superisoform_seq = None
    
    def __init__(self, ensg_id: str, biotype_filter=None):
        """
        Initialize a Gene object based on Ensembl ID

        params:

        """
        if biotype_filter is None:
            biotype_filter = ['protein_coding']
        self.ensg_id = ensg_id
        self.gene_info = self.download_gene_info()
        self.uniprot_id, self.refseq_id = self.check_alternate_id()
        self.start_pos, self.end_pos, self.strand = self.check_positions()
        self.transcripts = self.download_transcripts(self.ensg_id, biotype_filter)
        self.check_domain_redundancy()
        self.superisoform_seq = self.generate_superisofrom()

    def download_gene_info(self, ensg_id=None):
        if ensg_id is None:
            ensg_id = self.ensg_id
        mg = mygene.MyGeneInfo()
        get_gene_result = mg.getgene(ensg_id)
        return get_gene_result

    def download_transcripts(self, ensg_id=None, biotype_filter=None):
        if biotype_filter is None:
            biotype_filter = ['protein_coding']
        if ensg_id is None:
            ensg_id = self.ensg_id
        try:
            isoforms = ensembl_rest.lookup(ensg_id,
                                           params={'multiple_sequences': True,
                                                   'type': 'protein',
                                                   'expand': True
                                                   }
                                           )
        except ensembl_rest.HTTPError as err:
            error_code = err.response.status_code
            error_message= err.response.json()['error']
            if(error_code==400) and ("not found" in error_message):
                print(ensg_id + " not found")
            elif(error_code==400) and ("No sequences returned" in error_message):
                print(ensg_id + " no protein sequences found")
            else:
                raise
            return None
        transcripts = []
        for isoform in isoforms["Transcript"]:
            if isoform['biotype'] not in biotype_filter:
                continue
            transcript = Transcript(isoform['id'], isoform['Translation']['id'])
            transcripts.append(transcript)
        return transcripts

    def check_alternate_id(self, ensg_id = None):
        if ensg_id is None:
            ensg_id = self.ensg_id
        gene_info = self.gene_info
        print(gene_info['refseq'])
        uniprot_id = None
        refseq_id = None
        if gene_info is not None:
            if "uniprot" in gene_info:
                uniprot_id = gene_info["uniprot"]
            if "refseq" in gene_info:
                refseq_id = list(filter(re.compile("NG_.*").match, gene_info["refseq"]["genomic"][0]))
        if uniprot_id is None:
            print("No UniProt ID found for " + ensg_id)
        return uniprot_id, refseq_id

    def check_positions(self):
        gene_info = self.gene_info
        if type(gene_info['genomic_pos']) is list:
            start_pos = gene_info['genomic_pos'][0]['start']
            end_pos = gene_info['genomic_pos'][0]['end']
            strand = gene_info['genomic_pos'][0]['strand']
        else:
            start_pos = gene_info['genomic_pos']['start']
            end_pos = gene_info['genomic_pos']['end']
            strand = gene_info['genomic_pos']['strand']
        return start_pos, end_pos, strand

    def check_domain_redundancy(self, transcripts=None):
        if transcripts is None:
            transcripts = self.transcripts
        

    def generate_superisofrom(self):
        pass

