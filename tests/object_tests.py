import unittest
from Bio.SeqFeature import SimpleLocation

import bin.TF_ASIF.gene as gene
from bin.TF_ASIF.transcript import Transcript
from bin.TF_ASIF.domain import Domain


def _make_domain(types, start=0, end=10):
    """Build a Domain double with pre-set types, bypassing the real
    interpro_superfamily_domains_DBD.tsv-based classification lookup."""
    domain = Domain.__new__(Domain)
    domain.interpro_id = "IPR000000"
    domain.start = start
    domain.end = end
    domain.source = "SuperFamily"
    domain.pos = SimpleLocation(start, end)
    domain.types = types
    return domain


def _make_transcript(ensp_id, uniprot_id, domains):
    """Build a Transcript double with only the attributes
    check_domain_redundancy reads, bypassing the real constructor's
    file/network dependencies."""
    transcript = Transcript.__new__(Transcript)
    transcript.ensp_id = ensp_id
    transcript.uniprot_id = uniprot_id
    transcript.domains = domains
    return transcript


class TestCheckDomainRedundancy(unittest.TestCase):
    def test_domains_assigned_to_correct_transcript(self):
        # Non-overlapping coordinates so the (now default-on) overlap
        # merging can't collapse these into one another; this test is
        # specifically about prot_id attribution, not merging.
        domain_a = _make_domain(["DNA-binding"], start=0, end=10)
        domain_b = _make_domain(["DNA-binding"], start=20, end=30)
        domain_c = _make_domain(["DNA-binding"], start=40, end=50)
        transcript_a = _make_transcript("ENSP_A", "UNIPROT_A", [domain_a])
        transcript_b = _make_transcript("ENSP_B", "UNIPROT_B", [domain_b])
        transcript_c = _make_transcript("ENSP_C", "UNIPROT_C", [domain_c])

        test_gene = gene.Gene.__new__(gene.Gene)
        test_gene.transcripts = [transcript_a, transcript_b, transcript_c]
        test_gene.check_domain_redundancy()

        # Regression test: each transcript must get back only its own
        # domain(s), not another transcript's (see check_domain_redundancy
        # bugfix, 2026-09-24 - a stale loop variable used to misattribute
        # every pooled domain to whichever transcript was last in the list).
        self.assertEqual(transcript_a.filtered_domains, [domain_a])
        self.assertEqual(transcript_b.filtered_domains, [domain_b])
        self.assertEqual(transcript_c.filtered_domains, [domain_c])

    def test_transcript_without_uniprot_id_excluded(self):
        domain_a = _make_domain(["DNA-binding"])
        domain_b = _make_domain(["DNA-binding"])
        transcript_a = _make_transcript("ENSP_A", None, [domain_a])
        transcript_b = _make_transcript("ENSP_B", "UNIPROT_B", [domain_b])

        test_gene = gene.Gene.__new__(gene.Gene)
        test_gene.transcripts = [transcript_a, transcript_b]
        test_gene.check_domain_redundancy()

        self.assertEqual(transcript_a.filtered_domains, [])
        self.assertEqual(transcript_b.filtered_domains, [domain_b])

    def test_unclassified_domains_dropped(self):
        domain_dbd = _make_domain(["DNA-binding"])
        domain_other = _make_domain([])
        transcript_a = _make_transcript("ENSP_A", "UNIPROT_A", [domain_dbd, domain_other])

        test_gene = gene.Gene.__new__(gene.Gene)
        test_gene.transcripts = [transcript_a]
        test_gene.check_domain_redundancy()

        self.assertEqual(transcript_a.filtered_domains, [domain_dbd])

    def test_overlapping_domains_merged_by_default(self):
        domain_small = _make_domain(["DNA-binding"], start=0, end=10)
        domain_large = _make_domain(["DNA-binding"], start=5, end=20)
        transcript_small = _make_transcript("ENSP_SMALL", "UNIPROT_SMALL", [domain_small])
        transcript_large = _make_transcript("ENSP_LARGE", "UNIPROT_LARGE", [domain_large])

        test_gene = gene.Gene.__new__(gene.Gene)
        test_gene.transcripts = [transcript_small, transcript_large]
        test_gene.check_domain_redundancy()

        # merge_overlapping defaults to True: the two domains overlap (5
        # falls within 0-10), so only the larger one survives; the smaller
        # domain's transcript gets none.
        self.assertEqual(transcript_small.filtered_domains, [])
        self.assertEqual(transcript_large.filtered_domains, [domain_large])

    def test_keep_overlapping_domains_when_disabled(self):
        domain_small = _make_domain(["DNA-binding"], start=0, end=10)
        domain_large = _make_domain(["DNA-binding"], start=5, end=20)
        transcript_small = _make_transcript("ENSP_SMALL", "UNIPROT_SMALL", [domain_small])
        transcript_large = _make_transcript("ENSP_LARGE", "UNIPROT_LARGE", [domain_large])

        test_gene = gene.Gene.__new__(gene.Gene)
        test_gene.transcripts = [transcript_small, transcript_large]
        test_gene.check_domain_redundancy(merge_overlapping=False)

        # With merge_overlapping explicitly disabled, overlap between
        # transcripts is not collapsed; each transcript keeps its own domain.
        self.assertEqual(transcript_small.filtered_domains, [domain_small])
        self.assertEqual(transcript_large.filtered_domains, [domain_large])


class TestGene(unittest.TestCase):
    def test_Gene_creation(self):
        test_ensg_id = "ENSG00000101076"
        binding_site_file = "/mnt/data/storage/WPI/Korkin_Lab/DNA_Binding_IF/ppi_binding_sites.tsv"
        cds_fasta_file = "/mnt/data/storage/WPI/Korkin_Lab/DNA_Binding_IF/Homo_sapiens.GRCh38.cds.all.fa.gz"
        uniprot_mapping_file = "/mnt/data/storage/WPI/Korkin_Lab/DNA_Binding_IF/Homo_sapiens.GRCh38.109.uniprot.tsv.gz"
        gtf_file = "/mnt/data/storage/WPI/Korkin_Lab/DNA_Binding_IF/Homo_sapiens.GRCh38.109.gtf.gz"
        interpro_domains_file = "/mnt/data/storage/WPI/Korkin_Lab/DNA_Binding_IF/Homo_sapiens.GRCh38.interpro_domains.tsv.gz"
        test_gene = gene.Gene(test_ensg_id, binding_site_file, cds_fasta_file, uniprot_mapping_file,
                              gtf_file, interpro_domains_file)
        self.assertEqual(test_ensg_id, test_gene.ensg_id)  # add assertion here
        self.assertGreater(len(test_gene.transcripts), 0)
        self.assertIsNotNone(test_gene.symbol)
        self.assertGreater(test_gene.start_pos, 0)
        self.assertGreater(test_gene.end_pos, test_gene.start_pos)
        self.assertIn(test_gene.strand, [1, -1])
        self.assertEqual(test_gene.transcripts[0].enst_id[:4], "ENST")
        self.assertEqual(test_gene.transcripts[0].ensp_id[:4], "ENSP")
        self.assertGreater(len(test_gene.transcripts[0].domains), 0)
        test_enst_id = "ENST00000316673"
        expected_seq = ("MVSVNAPLGAPVESSYDTSPSEGTNLNAPNSLGVSALCAICGDRATGKHYGASSCDGCKG"
                        "FFRRSVRKNHMYSCRFSRQCVVDKDKRNQCRYCRLKKCFRAGMKKEAVQNERDRISTRRS"
                        "SYEDSSLPSINALLQAEVLSRQITSPVSGINGDIRAKKIASIADVCESMKEQLLVLVEWA"
                        "KYIPAFCELPLDDQVALLRAHAGEHLLLGATKRSMVFKDVLLLGNDYIVPRHCPELAEMS"
                        "RVSIRILDELVLPFQELQIDDNEYAYLKAIIFFDPDAKGLSDPGKIKRLRSQVQVSLEDY"
                        "INDRQYDSRGRFGELLLLLPTLQSITWQMIEQIQFIKLFGMAKIDNLLQEMLLGGSPSDA"
                        "PHAHHPLHPHLMQEHMGTNVIVANTMPTHLSNGQMCEWPRPRGQAATPETPQPSPPGGSG"
                        "SEPYKLLPGAVATIVKPLSAIPQPTITKQEVI")
        matching_transcripts = [t for t in test_gene.transcripts if t.enst_id == test_enst_id]
        self.assertEqual(len(matching_transcripts), 1)
        test_transcript = matching_transcripts[0]
        self.assertEqual(str(test_transcript.prot_seq).rstrip("*"), expected_seq)

        # GTF-derived CDS block structure (see Gene._get_gtf_cds_blocks): exon
        # lengths in transcript order, including the appended stop codon.
        expected_exons_rna = [SimpleLocation(0, 49), SimpleLocation(49, 224), SimpleLocation(224, 319),
                              SimpleLocation(319, 426), SimpleLocation(426, 582), SimpleLocation(582, 670),
                              SimpleLocation(670, 826), SimpleLocation(826, 1063), SimpleLocation(1063, 1216),
                              SimpleLocation(1216, 1356), SimpleLocation(1356, 1359)]
        expected_exons_prot = [SimpleLocation(0, 16), SimpleLocation(16, 74), SimpleLocation(74, 106),
                               SimpleLocation(106, 142), SimpleLocation(142, 194), SimpleLocation(194, 223),
                               SimpleLocation(223, 275), SimpleLocation(275, 354), SimpleLocation(354, 405),
                               SimpleLocation(405, 452), SimpleLocation(452, 453)]
        self.assertEqual(test_transcript.exons_rna, expected_exons_rna)
        self.assertEqual(test_transcript.exons_prot, expected_exons_prot)
        self.assertIsNotNone(test_transcript.rna_seq)
        self.assertEqual(len(test_transcript.rna_seq), test_transcript.exons_rna[-1].end)
        self.assertEqual(str(test_transcript.rna_seq.translate()), str(test_transcript.prot_seq))


if __name__ == '__main__':
    unittest.main()
