import unittest
import classes.gene as gene

class TestGene(unittest.TestCase):
    def test_Gene_creation(self):
        test_ensembl_id = "ENSG00000064489"
        test_gene = gene.Gene(test_ensembl_id)
        self.assertEqual(test_ensembl_id, test_gene.ensembl_id)  # add assertion here
        self.assertGreater(len(test_gene.download_transcripts()), 0)
        self.assertIsNotNone(test_gene.uniprot_id)

if __name__ == '__main__':
    unittest.main()
