import random

import biomart
import ensembl_rest
import mygene
import pandas as pd
from Bio import Entrez, SeqIO, SeqFeature, Align
from Bio.SeqFeature import SimpleLocation

Entrez.email = "smadha@wpi.edu"

# server = biomart.BiomartServer('http://may2025.archive.ensembl.org/biomart')
# mart = server.datasets['hsapiens_gene_ensembl']

#binding_site_df = pd.read_csv("../../../ppi_binding_sites.tsv", sep='\t', header=0)

class Domain:
    """
    Domain object
    """

    def __init__(self, interpro_id, start, end, source, pos=None, **kwargs):
        self.interpro_id = interpro_id
        self.start = start
        self.end = end
        self.source = source
        #self.types = self.determine_types()
        if pos is None:
            self.pos = SeqFeature.SeqFeature(SeqFeature.FeatureLocation(start, end))
        else:
            self.pos = pos

    def determine_types(self):
        types = []
        if self.determine_dna_binding():
            types.append("DNA-binding")
        if self.determine_protein_interaction():
            types.append("PPI")
        return types

    def determine_dna_binding(self, dna_binding_file="../../../interpro_superfamily_domains_DBD.tsv"):
        interpro_superfamily_domains_DBD = pd.read_csv(dna_binding_file, sep='\t', index_col=0)
        if self.interpro_id is None or self.source!="SuperFamily":
            return False
        return interpro_superfamily_domains_DBD.loc[self.interpro_id,"DNA-binding"]

    def determine_protein_interaction(self):
        return False

    def __repr__(self):
        return "Interpro ID %s at %s of types %s" % (self.interpro_id, self.pos, self.types)


class Transcript:
    """
    Transcript object
    """
    seq = None
    uniprot_id = None
    domains = []

    def __init__(self, gene, enst_id: str, ensp_id: str, domain_types):
        self.gene = gene
        self.enst_id = enst_id
        self.ensp_id = ensp_id
        # print("Ensembl ID: " + self.enst_id)
        self.uniprot_id = self.get_uniprot_id()
        # print("UniProt ID: " + str(self.uniprot_id))
        self.refseq_id = self.get_refseq_id()
        # print("RefSeq ID: " + str(self.refseq_id))
        self.seq = self.download_sequence()
        #print("Sequence: " + self.seq)
        #self.domains = self.download_domains()
        #print("Domains: " + str(self.domains))
        if 'ppi' in domain_types:
            self.domains = self.yue_ppi_locations()
        #print(len(self.domains))

    def download_sequence(self, ensp_id=None):
        if ensp_id==None:
            ensp_id = self.ensp_id
        seq = ensembl_rest.sequence_id(ensp_id)["seq"]
        return seq

    def get_uniprot_id(self):
        print(self.ensp_id)
        uniprot_ids = idmapping_df.loc[(idmapping_df[2].str.contains(self.ensp_id)) &
                                       (idmapping_df[1]=="Ensembl_PRO"), 0].tolist()
        if len(uniprot_ids)>0:
            print(uniprot_ids)
            return uniprot_ids[0]
        else:
            return None

    def get_refseq_id(self, uniprot_id = None):
        if uniprot_id is None:
            if self.uniprot_id is not None:
                uniprot_id = self.uniprot_id
            else:
                return None
        refseq_ids = idmapping_df.loc[(idmapping_df[0].str.contains(uniprot_id)) &
                                      (idmapping_df[1]=="RefSeq") &
                                      (idmapping_df[2].str.startswith("NP_")), 2].tolist()
        if len(refseq_ids)>0:
            print(refseq_ids[0])
            return refseq_ids[0]
        else:
            return None

    def download_domains(self, ensp_id=None):
        if ensp_id is None:
            ensp_id = self.ensp_id
        results = ensembl_rest.overlap_translation(ensp_id,
                                                   type="domain")
        return [Domain(interpro_id=res["interpro"], source = res["type"], **res) for res in results]

    def yue_ppi_locations(self):
        if self.uniprot_id is None:
            return
        if self.domains is None:
            domains = []
        else:
            domains = self.domains
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

    def align_to_reference(self, refmode="superisoform"):
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
            ref_seq = self.gene.transcripts[0].seq
        #print(ref_seq)
        transcript_seq = self.seq
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

class Gene:
    """
    Gene object
    """
    seq = None
    uniprot_id = None
    refseq_id = None
    transcripts = None
    superisoform_seq = None
    
    def __init__(self, ensg_id: str, binding_site_file, idmapping_file,
                 biotype_filter=None, refmode="superisoform", domain_filter=None):
        """
        Initialize a Gene object based on Ensembl ID

        params:

        """
        global binding_site_df
        binding_site_df = pd.read_csv(binding_site_file, sep='\t', header=0)
        global idmapping_df
        idmapping_df = pd.read_csv(idmapping_file, sep='\t', header=None)
        if biotype_filter is None:
            biotype_filter = ['protein_coding']
        if domain_filter is None:
            domain_filter = ['ppi']
        self.ensg_id = ensg_id
        self.gene_info = self.download_gene_info()
        self.uniprot_id, self.refseq_id, self.symbol = self.check_alternate_id()
        if self.refseq_id is None:
            return
        self.start_pos, self.end_pos, self.strand = self.check_positions()
        # print("downloading transcripts")
        self.transcripts = self.download_transcripts(self.ensg_id, biotype_filter=biotype_filter,
                                                     domain_types =domain_filter)
        # print("checking redundancy")
        self.check_domain_redundancy()
        # print("generating superisoform")
        if refmode == "superisoform":
            self.superisoform_seq, self.superdomains = self.generate_superisoform()
        # print(self.superisoform_seq)
        # print(self.superdomains)
        # print("aligning to superisoform")
        for transcript in self.transcripts:
            transcript.align_to_reference(refmode=refmode)

    def download_gene_info(self, ensg_id=None):
        if ensg_id is None:
            ensg_id = self.ensg_id

        # server = biomart.BiomartServer('http://useast.ensembl.org/biomart')
        # mart = server.datasets['hsapiens_gene_ensembl']
        #
        # attributes = ['ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id',
        #               'refseq_mrna', 'refseq_peptide']
        # response = mart.search({'attributes': attributes,
        #                         'filters': {'ensembl_gene_id': ensg_id}
        #                         })
        # data = response.raw.data.decode('ascii')
        #
        # id_df = pd.read_csv(StringIO(data), sep='\t')

        mg = mygene.MyGeneInfo()
        get_gene_result = mg.getgene(ensg_id)
        return get_gene_result

    def check_alternate_id(self, ensg_id = None):
        if ensg_id is None:
            ensg_id = self.ensg_id
        gene_info = self.gene_info
        uniprot_id = None
        refseq_id = None
        symbol = None
        if gene_info is not None:
            if "uniprot" in gene_info:
                uniprot_id = gene_info["uniprot"]
            if "refseq" in gene_info:
                refseq_id = gene_info["refseq"]["genomic"][0]
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
            i=0
            while i < len(gene_info['genomic_pos']):
                if gene_info['genomic_pos'][i]['ensemblgene']==self.ensg_id:
                    break
                i+=1
            if i==len(gene_info['genomic_pos']):
                i=0
            start_pos = gene_info['genomic_pos'][i]['start']
            end_pos = gene_info['genomic_pos'][i]['end']
            strand = gene_info['genomic_pos'][i]['strand']
        else:
            start_pos = gene_info['genomic_pos']['start']
            end_pos = gene_info['genomic_pos']['end']
            strand = gene_info['genomic_pos']['strand']
        return start_pos, end_pos, strand

    def download_transcripts(self, ensg_id=None, biotype_filter=None, domain_types=None):
        if biotype_filter is None:
            biotype_filter = ['protein_coding']
        if domain_types is None:
            domain_types = ['ppi']
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
            # try:
            transcript = Transcript(self, isoform['id'], isoform['Translation']['id'], domain_types)
            if transcript.refseq_id is not None and transcript.uniprot_id is not None:
                transcript.seq = transcript.download_sequence()
                #transcript.yue_ppi_locations()
                transcripts.append(transcript)
            # except:
            #     print("Error getting transcript: " + isoform['id'])
        return transcripts

    def check_domain_redundancy(self, transcripts=None):
        if transcripts is None:
            transcripts = self.transcripts
        keeping_domains = []
        for classification in ["DNA-binding", "PPI"]:
            domain_queue = []
            for transcript in transcripts:
                if transcript.refseq_id is None:
                    continue
                domain_queue += [domain for domain in transcript.domains
                                 if classification in domain.types]# and domain.source in ['SuperFamily', 'Yue']]
            # print("domain_queue:")
            # print(domain_queue)
            while len(domain_queue) > 0:
                currDomain = domain_queue.pop()
                # removeList = []
                # for i, domain in enumerate(domain_queue):
                #     if (currDomain.start <= domain.start <= currDomain.end) or \
                #             (currDomain.start <= domain.end <= currDomain.end) or \
                #             (currDomain.start >= domain.start and
                #              currDomain.end <= domain.end):
                #         if currDomain.end - currDomain.start < domain.end - domain.start:
                #             currDomain = domain
                #         removeList.append(i)
                # for i in removeList[-1::-1]:
                #     del domain_queue[i]
                currDomain.prot_id = transcript.refseq_id
                keeping_domains.append(currDomain)
        # print("keeping_domains:")
        # print(keeping_domains)
        for transcript in self.transcripts:
            transcript.filtered_domains = \
                [domain for domain in keeping_domains if domain.prot_id == transcript.refseq_id]

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
                refseq_id = feature.qualifiers['protein_id'][0]
                for transcript in self.transcripts:
                    # if transcript.refseq_id is None:
                    #     break
                    these_domains = []
                    if transcript.refseq_id in refseq_id:
                        these_domains = transcript.filtered_domains

                    parts = enumerate(feature.location.parts)
                    prepend_seq = ""
                    exon_start = 0
                    exon_end = 0
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
                                if type(domain.pos)==SimpleLocation:
                                    domain_start = None
                                    domain_end = None
                                    # domain in the exon
                                    if domain.start > exon_start and \
                                            domain.end < exon_end:
                                        domain_start = domain.start
                                        domain_end = domain.end
                                    # domain starts in the exon, but ends after
                                    elif exon_start < domain.start < exon_end < domain.end:
                                        domain_start = domain.start
                                        domain_end = exon_end
                                    # domain starts before exon, but ends in
                                    elif domain.start < exon_start < domain.end < exon_end:
                                        domain_start = exon_start
                                        domain_end = domain.end
                                    # domain contains exon
                                    elif domain.start < exon_start and \
                                            domain.end > exon_end:
                                        domain_start = exon_start
                                        domain_end = exon_end
                                    if domain_start is not None:
                                        exon_domains[".".join([str(domain_start),
                                                               str(domain_end),
                                                               str(max(domain.start - exon_start, 0)),
                                                               str(max(exon_end - domain.end, 0))])] = domain
                                else:
                                    domain_start = None
                                    domain_end = None
                                    # domain in the exon
                                    if domain.start > exon_start and \
                                            domain.end < exon_end:
                                        domain_start = domain.start
                                        domain_end = domain.end
                                    # domain starts in the exon, but ends after
                                    elif exon_start < domain.start < exon_end < domain.end:
                                        domain_start = domain.start
                                        domain_end = exon_end
                                    # domain starts before exon, but ends in
                                    elif domain.start < exon_start < domain.end < exon_end:
                                        domain_start = exon_start
                                        domain_end = domain.end
                                    # domain contains exon
                                    elif domain.start < exon_start and \
                                            domain.end > exon_end:
                                        domain_start = exon_start
                                        domain_end = exon_end
                                    if domain_start is not None:
                                        exon_domains[".".join([str(domain_start),
                                                               str(domain_end),
                                                               str(max(domain.start - exon_start, 0)),
                                                               str(max(exon_end - domain.end, 0))])] = domain
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
                                if type(domain.pos) == SimpleLocation:
                                    domain_start = None
                                    domain_end = None
                                    # domain in the exon
                                    if domain.start > exon_start and \
                                            domain.end < exon_end:
                                        domain_start = domain.start
                                        domain_end = domain.end
                                    # domain starts in the exon, but ends after
                                    elif exon_start < domain.start < exon_end < domain.end:
                                        domain_start = domain.start
                                        domain_end = exon_end
                                    # domain starts before exon, but ends in
                                    elif domain.start < exon_start < domain.end < exon_end:
                                        domain_start = exon_start
                                        domain_end = domain.end
                                    # domain contains exon
                                    elif domain.start < exon_start and \
                                            domain.end > exon_end:
                                        domain_start = exon_start
                                        domain_end = exon_end
                                    if domain_start is not None:
                                        exon_domains[".".join([str(domain_start).strip("<>"),
                                                               str(domain_end).strip("<>"),
                                                               str(max(domain.start - exon_start, 0)).strip("<>"),
                                                               str(max(exon_end - domain.end, 0)).strip("<>")])] = domain
                                else:
                                    domain_start = None
                                    domain_end = None
                                    # domain in the exon
                                    if domain.start > exon_start and \
                                            domain.end < exon_end:
                                        domain_start = domain.start
                                        domain_end = domain.end
                                    # domain starts in the exon, but ends after
                                    elif exon_start < domain.start < exon_end < domain.end:
                                        domain_start = domain.start
                                        domain_end = exon_end
                                    # domain starts before exon, but ends in
                                    elif domain.start < exon_start < domain.end < exon_end:
                                        domain_start = exon_start
                                        domain_end = domain.end
                                    # domain contains exon
                                    elif domain.start < exon_start and \
                                            domain.end > exon_end:
                                        domain_start = exon_start
                                        domain_end = exon_end
                                    if domain_start is not None:
                                        exon_domains[".".join([str(domain_start).strip("<>"),
                                                               str(domain_end).strip("<>"),
                                                               str(max(domain.start - exon_start, 0)).strip("<>"),
                                                               str(max(exon_end - domain.end, 0)).strip("<>")])] = domain

                        exon_key = ".".join([str(location.start).strip("<>"),
                                             str(location.end).strip("<>"),
                                             prepend_seq,
                                             str(tail_len)])
                        if not exon_key in reformed_exons:
                            reformed_exons[exon_key] = {'seq': this_seq.translate().seq, 'domains': exon_domains}
                        else:
                            reformed_exons[exon_key]['domains'].update(exon_domains)
                        prepend_seq = str(next_location.extract(record).seq)

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

        if strand == -1:
            tmp = tmp[::-1]
        superdomains = []
        for exon in tmp:
            if len(reformed_exons[exon]['domains']) > 0:
                for domain in reformed_exons[exon]['domains']:
                    domain_start = int(domain.split(".")[2]) + len(superisoform)
                    domain_end = int(domain.split(".")[1]) - int(domain.split(".")[0]) + domain_start
                    superdomains.append(Domain(interpro_id="SD" + str(random.randint(0, 9999)), source="Yue", start=domain_start,
                           end=domain_end, pos=SeqFeature.FeatureLocation(domain_start, domain_end+1)))
            superisoform += reformed_exons[exon]['seq']
        return superisoform, superdomains
