import gzip
import random
import re

import mygene
import pandas as pd
from Bio import SeqFeature
from Bio.SeqFeature import SimpleLocation

from bin.TF_ASIF.domain import Domain
from bin.TF_ASIF.transcript import Transcript

class Gene:
    """
    Gene object
    """
    seq = None
    uniprot_id = None
    refseq_id_chrom = None
    transcripts = []
    superisoform_seq = None
    _gtf_index = None
    _gtf_cds_blocks = None
    _gtf_exon_coords = None

    def __init__(self, ensg_id: str, binding_site_file, idmapping_file, cds_fasta_file, uniprot_mapping_file,
                 gtf_file, biotype_filter=None, refmode="superisoform", domain_filter=None):
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
            through to each Transcript for RefSeq ID resolution.
        cds_fasta_file : str
            Path to a (optionally gzipped) Ensembl "cds.all.fa" FASTA
            file, passed through to each Transcript for local protein
            sequence lookup.
        uniprot_mapping_file : str
            Path to a (optionally gzipped) Ensembl protein-to-UniProt
            xref TSV file (e.g. "*.uniprot.tsv.gz"), passed through to
            each Transcript for UniProt ID resolution.
        gtf_file : str
            Path to a (optionally gzipped) Ensembl GTF annotation
            file, used to list this gene's isoforms.
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
                                                     idmapping_df=idmapping_df, cds_fasta_file=cds_fasta_file,
                                                     uniprot_mapping_file=uniprot_mapping_file, gtf_file=gtf_file)
        # print("checking redundancy")
        self.check_domain_redundancy()
        # print("generating superisoform")
        if refmode == "superisoform":
            self.superisoform_seq, self.superdomains = self.generate_superisoform(gtf_file=gtf_file)
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
                             idmapping_df=None, cds_fasta_file=None, uniprot_mapping_file=None, gtf_file=None):
        """
        Look up this gene's isoforms in a local Ensembl GTF annotation
        file and build a Transcript object for each isoform that
        passes biotype_filter and has both a UniProt and RefSeq ID,
        filling in its RNA/protein sequence and exon locations from
        that same GTF file's CDS block structure and a local Ensembl
        CDS FASTA file. Isoforms whose protein sequence or CDS block
        structure cannot be determined (e.g. a partial CDS) are
        skipped.

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
            Transcript for RefSeq ID resolution.
        cds_fasta_file : str, optional
            Path to a (optionally gzipped) Ensembl "cds.all.fa" FASTA
            file, passed through to each Transcript for local protein
            sequence lookup.
        uniprot_mapping_file : str, optional
            Path to a (optionally gzipped) Ensembl protein-to-UniProt
            xref TSV file (e.g. "*.uniprot.tsv.gz"), passed through to
            each Transcript for UniProt ID resolution.
        gtf_file : str, optional
            Path to a (optionally gzipped) Ensembl GTF annotation
            file, used to list this gene's isoforms.

        Returns
        -------
        list[Transcript]
            Transcripts for this gene, or an empty list if the gene
            has no isoforms in gtf_file.
        """
        if biotype_filter is None:
            biotype_filter = ['protein_coding']
        if domain_types is None:
            domain_types = ['ppi']
        if ensg_id is None:
            ensg_id = self.ensg_id
        isoforms = self._get_gtf_index(gtf_file).get(ensg_id)
        if not isoforms:
            print(ensg_id + " not found")
            return []

        gtf_cds_blocks = self._get_gtf_cds_blocks(gtf_file)
        cds_sequences = Transcript._get_cds_sequences(cds_fasta_file)

        transcripts = []
        for isoform in isoforms:
            if isoform['biotype'] not in biotype_filter:
                continue
            transcript = Transcript(self, isoform['id'], isoform['protein_id'], domain_types=domain_types,
                                    binding_site_df=binding_site_df, idmapping_df=idmapping_df,
                                    cds_fasta_file=cds_fasta_file, uniprot_mapping_file=uniprot_mapping_file)
            if transcript.refseq_id is not None and transcript.uniprot_id is not None:
                exon_lengths = gtf_cds_blocks.get(transcript.enst_id)
                prot_seq = transcript.download_sequence()
                if exon_lengths is None or prot_seq is None:
                    continue

                exons_rna = []
                exons_prot = []
                exon_rna_start = 0
                for exon_len in exon_lengths:
                    exon_rna_end = exon_len + exon_rna_start
                    exons_rna += [SimpleLocation(exon_rna_start, exon_rna_end)]
                    exon_prot_start = (exon_rna_start - exon_rna_start % 3) // 3
                    exons_prot += [SimpleLocation(exon_prot_start, exon_rna_end // 3)]
                    exon_rna_start = exon_rna_end

                transcript.rna_seq = cds_sequences.get(transcript.enst_id)
                transcript.exons_rna = exons_rna
                transcript.exons_prot = exons_prot
                transcript.prot_seq = prot_seq
                transcripts.append(transcript)
        return transcripts

    @classmethod
    def _get_gtf_index(cls, gtf_file):
        """
        Load and cache, from an Ensembl GTF annotation file, each
        gene's transcripts, keyed by Ensembl Gene ID, along with each
        coding transcript's CDS block structure (see
        _get_gtf_cds_blocks and _get_gtf_exon_coords). The file is
        only parsed once per process; subsequent calls reuse the
        cache.

        Parameters
        ----------
        gtf_file: str
            Path to a (optionally gzipped) Ensembl GTF annotation
            file.

        Returns
        -------
        dict[str, list[dict]]
            For each Ensembl Gene ID, a list of {"id": <Ensembl
            Transcript ID>, "biotype": <transcript biotype>,
            "protein_id": <Ensembl Protein ID, or None if the
            transcript has no CDS>} dicts.
        """
        if cls._gtf_index is None:
            attribute_re = re.compile(r'(\w+) "([^"]*)"')
            transcript_info = {}
            transcript_protein_ids = {}
            cds_blocks = {}
            opener = gzip.open if gtf_file.endswith(".gz") else open
            with opener(gtf_file, "rt") as handle:
                for line in handle:
                    if line.startswith("#"):
                        continue
                    fields = line.rstrip("\n").split("\t")
                    feature = fields[2]
                    if feature == "transcript":
                        attributes = dict(attribute_re.findall(fields[8]))
                        transcript_info[attributes["transcript_id"]] = \
                            (attributes["gene_id"], attributes["transcript_biotype"])
                    elif feature == "CDS" or feature == "stop_codon":
                        attributes = dict(attribute_re.findall(fields[8]))
                        transcript_id = attributes["transcript_id"]
                        if feature == "CDS":
                            transcript_protein_ids.setdefault(transcript_id, attributes.get("protein_id"))
                        cds_blocks.setdefault(transcript_id, []).append(
                            (int(fields[3]), int(fields[4]), fields[6], fields[7]))
            gtf_index = {}
            for transcript_id, (gene_id, biotype) in transcript_info.items():
                gtf_index.setdefault(gene_id, []).append({
                    "id": transcript_id,
                    "biotype": biotype,
                    "protein_id": transcript_protein_ids.get(transcript_id),
                })
            gtf_cds_blocks = {}
            gtf_exon_coords = {}
            for transcript_id, blocks in cds_blocks.items():
                blocks = sorted(blocks, key=lambda block: block[0])
                if blocks[0][2] == "-":
                    blocks = blocks[::-1]
                lengths = [end - start + 1 for start, end, strand, frame in blocks]
                coords = [(start, end, frame) for start, end, strand, frame in blocks]
                # A transcript with an incomplete (5' truncated) CDS starts
                # mid-codon; Ensembl's CDS FASTA left-pads it with Ns to a
                # full codon, per the GTF frame of the first CDS block. This
                # padding has no genomic position of its own, so its coord
                # entry is None.
                first_frame = blocks[0][3]
                if first_frame.isdigit() and int(first_frame) != 0:
                    lengths.insert(0, 3 - int(first_frame))
                    coords.insert(0, None)
                gtf_cds_blocks[transcript_id] = lengths
                gtf_exon_coords[transcript_id] = coords
            cls._gtf_index = gtf_index
            cls._gtf_cds_blocks = gtf_cds_blocks
            cls._gtf_exon_coords = gtf_exon_coords
        return cls._gtf_index

    @classmethod
    def _get_gtf_cds_blocks(cls, gtf_file):
        """
        Load and cache, from an Ensembl GTF annotation file, each
        coding transcript's CDS block lengths (including its stop
        codon), ordered in transcript (5' to 3') order. Concatenating
        a transcript's CDS blocks in this order and translating them
        reproduces the sequence in a matching Ensembl "cds.all.fa"
        FASTA file. The file is only parsed once per process (shared
        with _get_gtf_index); subsequent calls reuse the cache.

        Parameters
        ----------
        gtf_file: str
            Path to a (optionally gzipped) Ensembl GTF annotation
            file.

        Returns
        -------
        dict[str, list[int]]
            CDS block nucleotide lengths, in transcript order, keyed
            by Ensembl Transcript ID.
        """
        cls._get_gtf_index(gtf_file)
        return cls._gtf_cds_blocks

    @classmethod
    def _get_gtf_exon_coords(cls, gtf_file):
        """
        Load and cache, from an Ensembl GTF annotation file, each
        coding transcript's CDS block genomic coordinates and reading
        frame, index-aligned with _get_gtf_cds_blocks (and so also
        with the exons_rna/exons_prot built from it). The synthetic
        5' N-padding block inserted for an incomplete CDS (see
        _get_gtf_cds_blocks) has no genomic position and is
        represented as None. The file is only parsed once per process
        (shared with _get_gtf_index); subsequent calls reuse the
        cache.

        Parameters
        ----------
        gtf_file: str
            Path to a (optionally gzipped) Ensembl GTF annotation
            file.

        Returns
        -------
        dict[str, list[tuple | None]]
            For each Ensembl Transcript ID, a list of (genomic start,
            genomic end, reading frame) tuples (or None for the
            synthetic N-padding block), in transcript order.
        """
        cls._get_gtf_index(gtf_file)
        return cls._gtf_exon_coords

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

    def generate_superisoform(self, gtf_file=None):
        """
        Build a superisoform protein sequence by combining each unique
        exon across all of this gene's transcripts, deduplicated by
        genomic position and reading frame (so two transcripts sharing
        the same physical, in-frame exon contribute it only once) and
        arranged in genomic order, and remap each transcript's
        filtered domains onto the resulting superisoform coordinates.

        Parameters
        ----------
        gtf_file : str, optional
            Path to a (optionally gzipped) Ensembl GTF annotation
            file, used to look up each transcript's exon genomic
            coordinates and reading frame for deduplication and
            ordering.

        Returns
        -------
        tuple
            (superisoform, superdomains) where superisoform is the
            combined-exon protein sequence (str) and superdomains is
            the list[Domain] of domains mapped onto it.
        """
        exon_coords = self._get_gtf_exon_coords(gtf_file)

        unique_exons = {}
        for transcript in self.transcripts:
            if transcript.refseq_id is None or transcript.prot_seq is None:
                continue
            coords = exon_coords.get(transcript.enst_id)
            if coords is None:
                continue
            for exon, coord in zip(transcript.exons_prot, coords):
                # The synthetic N-padding block for an incomplete CDS (see
                # _get_gtf_cds_blocks) has no genomic position, so it can't
                # be deduplicated against other transcripts; skip it.
                if coord is None:
                    continue
                if coord not in unique_exons:
                    unique_exons[coord] = {"seq": str(exon.extract(transcript.prot_seq)), "domains": {}}

                for domain in transcript.filtered_domains:
                    domain_start = None
                    domain_end = None
                    # domain in the exon
                    if domain.start > exon.start and domain.end < exon.end:
                        domain_start = domain.start
                        domain_end = domain.end
                    # domain starts in the exon, but ends after
                    elif exon.start < domain.start < exon.end < domain.end:
                        domain_start = domain.start
                        domain_end = exon.end
                    # domain starts before exon, but ends in
                    elif domain.start < exon.start < domain.end < exon.end:
                        domain_start = exon.start
                        domain_end = domain.end
                    # domain contains exon
                    elif domain.start < exon.start and domain.end > exon.end:
                        domain_start = exon.start
                        domain_end = exon.end
                    if domain_start is not None:
                        domain_key = (domain_start, domain_end,
                                     max(domain.start - exon.start, 0),
                                     max(exon.end - domain.end, 0))
                        unique_exons[coord]["domains"][domain_key] = domain

        ordered_coords = sorted(unique_exons.keys(), key=lambda coord: coord[0])
        if self.strand == -1:
            ordered_coords = ordered_coords[::-1]

        superisoform = ""
        superdomains = []
        for coord in ordered_coords:
            exon = unique_exons[coord]
            for domain_start, domain_end, offset_before, offset_after in exon["domains"]:
                remapped_start = offset_before + len(superisoform)
                remapped_end = domain_end - domain_start + remapped_start
                superdomains.append(Domain(interpro_id="SD" + str(random.randint(0, 9999)), source="Yue",
                                           start=remapped_start, end=remapped_end,
                                           pos=SeqFeature.FeatureLocation(remapped_start, remapped_end + 1)))
            superisoform += exon["seq"]
        return superisoform, superdomains
