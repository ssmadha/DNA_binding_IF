import random

import ensembl_rest
import mygene
import pandas as pd
from Bio import Entrez, SeqIO, SeqFeature
from Bio.SeqFeature import SimpleLocation

from bin.TF_ASIF.domain import Domain
from bin.TF_ASIF.transcript import Transcript

Entrez.email = "smadha@wpi.edu"

class Gene:
    """
    Gene object
    """
    seq = None
    uniprot_id = None
    refseq_id_chrom = None
    transcripts = []
    superisoform_seq = None
    
    def __init__(self, ensg_id: str, binding_site_file, idmapping_file,
                 biotype_filter=None, refmode="superisoform", domain_filter=None):
        """
        Constructor

        Downloads gene and transcript information, builds the transcripts'
        domains, and (if requested) generates and aligns a superisoform
        reference.

        Parameters
        ----------
        ensg_id : str
            Ensembl Gene ID.
        binding_site_file : str
            Path to a tab-separated file of PPI binding sites, passed
            through to each Transcript.
        idmapping_file : str
            Path to a tab-separated UniProt ID mapping file, passed
            through to each Transcript.
        biotype_filter : list, optional
            Ensembl transcript biotypes to keep. Defaults to
            ['protein_coding'].
        refmode : str, optional
            Reference mode used to align transcripts. "superisoform"
            generates a combined-exon reference sequence; any other
            value aligns to the first transcript instead.
        domain_filter : list, optional
            Domain types to include for each transcript. Defaults to
            ["ppi_domain", "dbi"].
        """
        binding_site_df = pd.read_csv(binding_site_file, sep='\t', header=0)
        idmapping_df = pd.read_csv(idmapping_file, sep='\t', header=None)
        if biotype_filter is None:
            biotype_filter = ['protein_coding']
        if domain_filter is None:
            domain_filter = ["ppi_domain", "dbi"]
        self.ensg_id = ensg_id
        self.gene_info = self.download_gene_info()
        self.uniprot_id, self.refseq_id_chrom, self.symbol = self.check_alternate_id()
        if self.refseq_id_chrom is None:
            return
        self.start_pos, self.end_pos, self.strand = self.check_positions()
        # print("downloading transcripts")
        self.transcripts = self.download_transcripts(self.ensg_id, biotype_filter=biotype_filter,
                                                     domain_types=domain_filter, binding_site_df=binding_site_df,
                                                     idmapping_df=idmapping_df)
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
        """
        Download gene information from MyGene.info

        Parameters
        ----------
        ensg_id : str, optional
            Ensembl Gene ID. Defaults to self.ensg_id.

        Returns
        -------
        dict
            Gene info result returned by mygene.MyGeneInfo.getgene.
        """
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
        """
        Identify the UniProt ID, RefSeq chromosome ID, and gene symbol
        for this gene from its downloaded gene info.

        Parameters
        ----------
        ensg_id : str, optional
            Ensembl Gene ID, used only for logging when an ID is
            missing. Defaults to self.ensg_id.

        Returns
        -------
        tuple
            (uniprot_id, refseq_id_chrom, symbol), any of which may be
            None if not found in gene_info.
        """
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
        """
        Determine the genomic start, end, and strand of this gene from
        gene_info, selecting the entry matching self.ensg_id when
        genomic_pos contains multiple entries.

        Returns
        -------
        tuple
            (start_pos, end_pos, strand).
        """
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

    def download_transcripts(self, ensg_id=None, biotype_filter=None, domain_types=None, binding_site_df=None,
                             idmapping_df=None):
        """
        Look up this gene's isoforms via the Ensembl REST API and build
        a Transcript object for each isoform that passes biotype_filter
        and has both a UniProt and RefSeq ID, filling in its RNA/protein
        sequence and exon locations from the matching RefSeq CDS
        feature (or by downloading the sequence if no match is found).

        Parameters
        ----------
        ensg_id : str, optional
            Ensembl Gene ID. Defaults to self.ensg_id.
        biotype_filter : list, optional
            Ensembl transcript biotypes to keep. Defaults to
            ['protein_coding'].
        domain_types : list, optional
            Domain types to pass through to each Transcript. Defaults
            to ['ppi'].
        binding_site_df : pd.DataFrame, optional
            Data frame of PPI binding sites, passed through to each
            Transcript.
        idmapping_df : pd.DataFrame, optional
            UniProt ID mapping data frame, passed through to each
            Transcript.

        Returns
        -------
        list[Transcript]
            Transcripts for this gene, or an empty list if the gene or
            its protein sequences could not be found on Ensembl.
        """
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
            return []

        refseq_id_chrom, start_pos, end_pos, strand = self.refseq_id_chrom, self.start_pos, self.end_pos, self.strand
        handle = Entrez.efetch(db="nucleotide",
                               id=refseq_id_chrom,
                               seq_start=start_pos,
                               seq_stop=end_pos,
                               rettype="gb")
        record = SeqIO.read(handle, "gb")
        handle.close()

        features = [feature for feature in record.features if feature.type=="CDS"]

        transcripts = []
        for isoform in isoforms["Transcript"]:
            if isoform['biotype'] not in biotype_filter:
                continue
            # try:
            transcript = Transcript(self, isoform['id'], isoform['Translation']['id'], domain_types=domain_types,
                                    binding_site_df=binding_site_df, idmapping_df=idmapping_df)
            if transcript.refseq_id is not None and transcript.uniprot_id is not None:
                prot_seq = None
                for feature in features:
                    if transcript.refseq_id in feature.qualifiers['protein_id'][0]:
                        rna_seq = feature.location.extract(record).seq
                        prot_seq = rna_seq.translate()

                        exons_rna = []
                        exons_prot = []
                        exon_rna_start = 0
                        for part in feature.location.parts:
                            exon_rna_end = len(part) + exon_rna_start
                            exons_rna += [SimpleLocation(exon_rna_start, exon_rna_end)]
                            exon_prot_start = (exon_rna_start - exon_rna_start%3)//3
                            exons_prot += [SimpleLocation(exon_prot_start, exon_rna_end//3)]
                            exon_rna_start = exon_rna_end

                        transcript.rna_seq = rna_seq
                        transcript.exons_rna = exons_rna
                        transcript.exons_prot = exons_prot
                        break

                if prot_seq is None:
                    prot_seq = transcript.download_sequence()
                transcript.prot_seq = prot_seq
                #transcript.yue_ppi_locations()
                transcripts.append(transcript)
            # except Exception as err:
            #     print("Error getting transcript: " + isoform['id'], "with error", err, file=sys.stderr)
        return transcripts

    def check_domain_redundancy(self, transcripts=None):
        """
        Collapse redundant domains across transcripts by classification
        (DNA-binding, PPI), then assign the surviving domains back to
        each transcript as transcript.filtered_domains.

        Parameters
        ----------
        transcripts : list[Transcript], optional
            Transcripts to check for domain redundancy. Defaults to
            self.transcripts.
        """
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
        """
        Build a superisoform protein sequence by combining each unique
        exon across all of this gene's transcripts, in genomic order,
        and remap each transcript's filtered domains onto the
        resulting superisoform coordinates.

        Returns
        -------
        tuple
            (superisoform, superdomains) where superisoform is the
            combined-exon protein sequence (str) and superdomains is
            the list[Domain] of domains mapped onto it.
        """
        superisoform_exon = ""
        superisoform_exons = []
        for transcript in self.transcripts:
            if transcript.refseq_id is not None:
                for exon in transcript.exons_prot:
                    exon_seq = str(exon.extract(transcript.prot_seq))
                    if exon_seq not in superisoform_exon:
                        print(transcript.domains)
                        superisoform_exon += exon_seq
                        if len(superisoform_exons)==0:
                            superisoform_exons.append(exon)
                        else:
                            superisoform_exons.append(SimpleLocation(superisoform_exons[-1].end,
                                                                     superisoform_exons[-1].end + exon.end - exon.start))
                        print(superisoform_exons)


        refseq_id_chrom = self.refseq_id_chrom
        symbol = self.symbol
        start_pos = self.start_pos
        end_pos = self.end_pos
        strand = self.strand
        handle = Entrez.efetch(db="nucleotide",
                               id=refseq_id_chrom,
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
