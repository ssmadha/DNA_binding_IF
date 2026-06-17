import unittest
import bin.TF_ASIF.gene as gene

class TestGene(unittest.TestCase):
    def test_Gene_creation(self):
        test_ensg_id = "ENSG00000068323"
        test_gene = gene.Gene(test_ensg_id)
        self.assertEqual(test_ensg_id, test_gene.ensg_id)  # add assertion here
        self.assertGreater(len(test_gene.download_transcripts()), 0)
        self.assertIsNotNone(test_gene.uniprot_id)
        self.assertIsNotNone(test_gene.refseq_id_chrom)
        self.assertIsNotNone(test_gene.symbol)
        self.assertGreater(test_gene.start_pos, 0)
        self.assertGreater(test_gene.end_pos, test_gene.start_pos)
        self.assertIn(test_gene.strand, [1, -1])
        self.assertEqual(test_gene.transcripts[0].enst_id[:4], "ENST")
        self.assertEqual(test_gene.transcripts[0].ensp_id[:4], "ENSP")
        self.assertGreater(len(test_gene.transcripts[0].domains), 0)


if __name__ == '__main__':
    unittest.main()
