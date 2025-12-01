import unittest
import classes.gene as gene

class TestGene(unittest.TestCase):
    test_ensg_id = "ENSG00000004848"
    test_gene = gene.Gene(test_ensg_id)
    def test_Gene_creation(self):
        self.assertEqual(self.test_ensg_id, self.test_gene.ensg_id)  # add assertion here
        self.assertGreater(len(self.test_gene.download_transcripts()), 0)
        self.assertIsNotNone(self.test_gene.uniprot_id)
        self.assertIsNotNone(self.test_gene.refseq_id)
        self.assertGreater(self.test_gene.start_pos, 0)
        self.assertGreater(self.test_gene.end_pos, self.test_gene.start_pos)
        self.assertIn(self.test_gene.strand, [1, -1])

    def test_Transcript_creation(self):
        self.assertEqual(self.test_gene.transcripts[0].enst_id[:4], "ENST")
        self.assertEqual(self.test_gene.transcripts[0].ensp_id[:4], "ENSP")
        self.assertGreater(len(self.test_gene.transcripts[0].domains), 0)
        self.assertIsNotNone(self.test_gene.transcripts[0].domains[0].interpro_id)

if __name__ == '__main__':
    unittest.main()
