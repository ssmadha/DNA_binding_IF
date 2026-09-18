import unittest
from Bio.SeqFeature import SimpleLocation

import bin.TF_ASIF.gene as gene

class TestGene(unittest.TestCase):
    def test_Gene_creation(self):
        test_ensg_id = "ENSG00000101076"
        binding_site_file = "/mnt/data/storage/WPI/Korkin_Lab/DNA_Binding_IF/ppi_binding_sites.tsv"
        idmapping_file = "/mnt/data/storage/WPI/Korkin_Lab/DNA_Binding_IF/HUMAN_9606_idmapping.filtered.dat"
        cds_fasta_file = "/mnt/data/storage/WPI/Korkin_Lab/DNA_Binding_IF/Homo_sapiens.GRCh38.cds.all.fa.gz"
        uniprot_mapping_file = "/mnt/data/storage/WPI/Korkin_Lab/DNA_Binding_IF/Homo_sapiens.GRCh38.109.uniprot.tsv.gz"
        gtf_file = "/mnt/data/storage/WPI/Korkin_Lab/DNA_Binding_IF/Homo_sapiens.GRCh38.109.gtf.gz"
        test_gene = gene.Gene(test_ensg_id, binding_site_file, idmapping_file, cds_fasta_file, uniprot_mapping_file,
                              gtf_file)
        self.assertEqual(test_ensg_id, test_gene.ensg_id)  # add assertion here
        self.assertGreater(len(test_gene.transcripts), 0)
        self.assertIsNotNone(test_gene.uniprot_id)
        self.assertIsNotNone(test_gene.refseq_id_chrom)
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
