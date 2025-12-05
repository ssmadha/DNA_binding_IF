
import ensembl_rest
import mygene
import re

import pandas as pd
from Bio import Entrez, SeqIO


class Domain:
    """
    Domain object
    """

    def __init__(self, interpro_id, start, end, source, **kwargs):
        self.interpro_id = interpro_id
        self.start = start
        self.end = end
        self.source = source
        self.types = self.determine_types()

    def determine_types(self):
        types = []
        if self.determine_dna_binding():
            types.append("DNA-binding")
        if self.determine_protein_interaction():
            types.append("Protein-interaction")

    def determine_dna_binding(self, dna_binding_file="interpro_superfamily_domains_DBD.tsv"):
        interpro_superfamily_domains_DBD = pd.read_csv("interpro_superfamily_domains_DBD.tsv", sep='\t', index_col=0)
        if self.interpro_id is None or self.source!="SuperFamily":
            return False
        return interpro_superfamily_domains_DBD.loc[self.interpro_id,"DNA-binding"]

    def determine_protein_interaction(self):
        pass


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
        self.uniprot_id, self.refseq_id, self.symbol = self.check_alternate_id()
        self.start_pos, self.end_pos, self.strand = self.check_positions()
        self.transcripts = self.download_transcripts(self.ensg_id, biotype_filter)
        self.check_domain_redundancy()
        self.superisoform_seq = self.generate_superisoform()

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
        symbol = None
        if gene_info is not None:
            if "uniprot" in gene_info:
                uniprot_id = gene_info["uniprot"]
            if "refseq" in gene_info:
                refseq_id = list(filter(re.compile("NG_.*").match, gene_info["refseq"]["genomic"][0]))
            if "symbol" in gene_info:
                symbol = gene_info["symbol"]
        if uniprot_id is None:
            print("No UniProt ID found for " + ensg_id)
        if refseq_id is None:
            print("No RefSeq ID found for " + ensg_id)
        if symbol is None:
            print("No gene symbol found for " + ensg_id)
        return uniprot_id, refseq_id, symbol

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
        

    def generate_superisoform(self):
        refseq_id = self.refseq_id
        symbol = self.symbol
        start_pos = self.start_pos
        end_pos = self.end_pos
        strand = self.strand
        handle = Entrez.efetch(db="nucleotide",
                               id=refseq_id,
                               seq_start=start_pos,
                               seq_stop=end_pos,
                               rettype="gb")
        record = SeqIO.read(handle, "gb")
        handle.close()
        reformed_exons = {}

        for feature in record.features:
            if feature.type == "CDS" and symbol in feature.qualifiers['gene']:
                # print(feature)
                refseq_id = feature.qualifiers['protein_id'][0]
                for transcript in self.transcripts:
                    if "refseq_id" in DNAbinding_seqs[gene_id][prot_id] and \
                            DNAbinding_seqs[gene_id][prot_id]['refseq_id'] != "" and \
                            re.search(DNAbinding_seqs[gene_id][prot_id]['refseq_id'], refseq_id) is not None:
                        break
                these_domains = []
                if 'refseq_id' in DNAbinding_seqs[gene_id][prot_id] and \
                        DNAbinding_seqs[gene_id][prot_id]['refseq_id'] != "" and \
                        re.search(DNAbinding_seqs[gene_id][prot_id]['refseq_id'], refseq_id):
                    these_domains = DNAbinding_seqs[gene_id][prot_id]['filtered_domains']
                parts = enumerate(feature.location.parts)
                prepend_seq = ""
                exon_start = 0
                exon_end = 0
                # print(prot_id)
                # print(these_domains)
                for i, location in parts:
                    exon_domains = {}
                    tail_len = (len(prepend_seq) + len(location)) % 3
                    if location.strand == 1:
                        this_location = SeqFeature.FeatureLocation(
                            location.start,
                            location.end - tail_len,
                            strand=location.strand)
                        next_location = SeqFeature.FeatureLocation(
                            location.end - tail_len, location.end,
                            strand=location.strand)
                        this_seq = prepend_seq + this_location.extract(record)
                        exon_len = int(len(this_seq) / 3)
                        exon_start = exon_end
                        exon_end = exon_start + exon_len
                        for domain in these_domains:
                            domain_start = None
                            domain_end = None
                            # domain in the exon
                            if domain['start'] > exon_start and \
                                    domain['end'] < exon_end:
                                domain_start = domain['start']
                                domain_end = domain['end']
                            # domain starts in the exon, but ends after
                            elif domain['start'] > exon_start and \
                                    domain['start'] < exon_end and \
                                    domain['end'] > exon_end:
                                domain_start = domain['start']
                                domain_end = exon_end
                            # domain starts before exon, but ends in
                            elif domain['start'] < exon_start and \
                                    domain['end'] > exon_start and \
                                    domain['end'] < exon_end:
                                domain_start = exon_start
                                domain_end = domain['end']
                            # domain contains exon
                            elif domain['start'] < exon_start and \
                                    domain['end'] > exon_end:
                                domain_start = exon_start
                                domain_end = exon_end
                            if domain_start is not None:
                                # print(domain)
                                # print(domain['start'], exon_start, domain_start, domain['start'] - exon_start)
                                # print(domain['end'], exon_end, domain_end, exon_end - domain['end'])
                                exon_domains[".".join([str(domain_start),
                                                       str(domain_end),
                                                       str(max(domain['start'] - exon_start, 0)),
                                                       str(max(exon_end - domain['end'], 0))])] = domain
                    elif location.strand == -1:
                        this_location = SeqFeature.FeatureLocation(
                            location.start + tail_len,
                            location.end,
                            strand=location.strand)
                        next_location = SeqFeature.FeatureLocation(
                            location.start, location.start + tail_len,
                            strand=location.strand)
                        this_seq = prepend_seq + this_location.extract(record)
                        exon_len = int(len(this_seq) / 3)
                        exon_start = exon_end
                        exon_end = exon_start + exon_len
                        for domain in these_domains:
                            domain_start = None
                            domain_end = None
                            # domain in the exon
                            if domain['start'] > exon_start and \
                                    domain['end'] < exon_end:
                                domain_start = domain['start']
                                domain_end = domain['end']
                            # domain starts in the exon, but ends after
                            elif domain['start'] > exon_start and \
                                    domain['start'] < exon_end and \
                                    domain['end'] > exon_end:
                                domain_start = domain['start']
                                domain_end = exon_end
                            # domain starts before exon, but ends in
                            elif domain['start'] < exon_start and \
                                    domain['end'] > exon_start and \
                                    domain['end'] < exon_end:
                                domain_start = exon_start
                                domain_end = domain['end']
                            # domain contains exon
                            elif domain['start'] < exon_start and \
                                    domain['end'] > exon_end:
                                domain_start = exon_start
                                domain_end = exon_end
                            if domain_start is not None:
                                # print(domain)
                                # print(domain['start'], exon_start, domain_start, domain['start'] - exon_start)
                                # print(domain['end'], exon_end, domain_end, exon_end - domain['end'])
                                exon_domains[".".join([str(domain_start).strip("<>"),
                                                       str(domain_end).strip("<>"),
                                                       str(max(domain['start'] - exon_start, 0)).strip("<>"),
                                                       str(max(exon_end - domain['end'], 0)).strip("<>")])] = domain
                    exon_key = ".".join([str(location.start).strip("<>"),
                                         str(location.end).strip("<>"),
                                         prepend_seq,
                                         str(tail_len)])
                    if not exon_key in reformed_exons:
                        reformed_exons[exon_key] = {'seq': this_seq.translate().seq, 'domains': exon_domains}
                    else:
                        reformed_exons[exon_key]['domains'].update(exon_domains)
                    prepend_seq = str(next_location.extract(record).seq)
                    # if isinstance(location.start, SeqFeature.BeforePosition):
                    # print(location.start)
        # print(reformed_exons)
        tmp = [el.split(".") for el in list(reformed_exons.keys())]
        for el in tmp:
            el[0] = int(el[0])
            el[1] = int(el[1])
            el[3] = int(el[3])
        tmp.sort()
        for el in tmp:
            el[0] = str(el[0])
            el[1] = str(el[1])
            el[3] = str(el[3])
        tmp = [".".join(el) for el in tmp]
        superisoform = ""
        # print(tmp)
        if strand == -1:
            tmp = tmp[::-1]
        superdomains = {}
        for exon in tmp:
            if len(reformed_exons[exon]['domains']) > 0:
                # print(reformed_exons[exon]['domains'])
                for domain in reformed_exons[exon]['domains']:
                    domain_start = int(domain.split(".")[2]) + len(superisoform)
                    domain_end = int(domain.split(".")[1]) - int(domain.split(".")[0]) + domain_start
                    # print(reformed_exons[exon]['domains'][domain]['seq_region_name'], domain_start, domain_end)
                    if str(reformed_exons[exon]['domains'][domain]) not in superdomains:
                        superdomains[str(reformed_exons[exon]['domains'][domain]['interpro'])] = \
                            SeqFeature.FeatureLocation(domain_start, domain_end)
                    else:
                        superdomains[str(reformed_exons[exon]['domains'][domain]['interpro'])] += \
                            SeqFeature.FeatureLocation(domain_start, domain_end)
            superisoform += reformed_exons[exon]['seq']
            # print(exon)
            # print(reformed_exons[exon]['seq'])
        DNAbinding_seqs[gene_id]['superisoform'] = superisoform
        DNAbinding_seqs[gene_id]['superdomains'] = superdomains
