import unittest
import warnings
from Bio.SeqFeature import SimpleLocation

import bin.TF_ASIF.gene as gene
from bin.TF_ASIF.transcript import Transcript
from bin.TF_ASIF.domain import Domain
from bin.TF_ASIF.superisoform import Superisoform
from bin.TF_ASIF import matching
from Bio.Seq import Seq


def _make_domain(types, start=0, end=10):
    """Build a Domain double with pre-set types, bypassing the real
    interpro_superfamily_domains_DBD.tsv-based classification lookup."""
    domain = Domain.__new__(Domain)
    domain.domain_id = "IPR000000"
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


class TestAlignToReference(unittest.TestCase):
    def _coverage(self, transcript_seq, identical_only, reference="ACDEFGHIK", location=SimpleLocation(0, 9)):
        # One superdomain, by default over the whole 9-residue reference.
        test_gene = gene.Gene.__new__(gene.Gene)
        test_gene.ensg_id = "G"
        test_gene.superisoform_seq = reference
        test_gene.superdomains = [Domain("SD1", location, "Yue")]
        transcript = Transcript.__new__(Transcript)
        transcript.gene = test_gene
        transcript.enst_id = "T"
        transcript.prot_seq = transcript_seq
        return transcript.align_to_reference(identical_only=identical_only)["SD1"]

    def test_mismatch_counts_as_covered_by_default(self):
        self.assertEqual(self._coverage("ACDWFGHIK", identical_only=False), 1.0)

    def test_insertion_does_not_shift_domain_positions(self):
        # The transcript has 3 extra residues after "C" and lacks the
        # reference's last 3 ("LMN"), so a domain over LMN is uncovered.
        # Indexing by alignment column would instead read "HIK" there.
        coverage = self._coverage("ACWWWDEFGHIK", identical_only=False,
                                  reference="ACDEFGHIKLMN", location=SimpleLocation(9, 12))
        self.assertEqual(coverage, 0.0)

    def test_mismatch_uncovered_when_identical_only(self):
        self.assertAlmostEqual(self._coverage("ACDWFGHIK", identical_only=True), 8 / 9)


class TestSuperisoform(unittest.TestCase):
    def test_alt_splice_site_shares_segments(self):
        # T2's second exon starts 6 nt later (in frame), so the shared
        # 206-220 part must be a single segment used by both.
        si = Superisoform("G", 1, {"T1": [(100, 109, "0"), (200, 220, "2")],
                                   "T2": [(100, 109, "0"), (206, 220, "2")]})
        self.assertEqual([s.key for s in si.segments],
                         [(100, 109, 0), (200, 205, 2), (206, 220, 2)])
        self.assertEqual(si.transcript_segments, {"T1": [0, 1, 2], "T2": [0, 2]})
        self.assertEqual(si.segments[1].transcripts, {"T1"})
        self.assertEqual(si.segments[2].transcripts, {"T1", "T2"})

    def test_same_span_different_phase_is_different_segment(self):
        # T2's first exon is 1 nt shorter, shifting the frame of the
        # shared 200-220 span.
        si = Superisoform("G", 1, {"T1": [(100, 109, "0"), (200, 220, "2")],
                                   "T2": [(100, 108, "0"), (200, 220, "0")]})
        self.assertEqual([s.key for s in si.segments],
                         [(100, 108, 0), (109, 109, 0), (200, 220, 0), (200, 220, 2)])
        self.assertEqual(si.transcript_segments, {"T1": [0, 1, 3], "T2": [0, 2]})

    def test_minus_strand_order(self):
        # On the minus strand a block's 5' end is its highest coordinate:
        # T2 skips 112-120 (9 nt, in frame), so 100-111 is shared.
        si = Superisoform("G", -1, {"T1": [(300, 311, "0"), (100, 120, "0")],
                                    "T2": [(300, 311, "0"), (100, 111, "0")]})
        self.assertEqual([s.key for s in si.segments],
                         [(300, 311, 0), (112, 120, 0), (100, 111, 0)])
        self.assertEqual(si.transcript_segments, {"T1": [0, 1, 2], "T2": [0, 2]})

    def test_residue_segments_split_codon(self):
        si = Superisoform("G", 1, {"T1": [(100, 109, "0"), (200, 205, "2")]})
        # 16 nt -> 5 complete codons; codon 4 (nt 9-11) spans both segments
        self.assertEqual(si.residue_segments("T1"), [(0,), (0,), (0,), (0, 1), (1,)])

    def test_incomplete_5_prime_padding(self):
        # GTF frame 1 -> Ensembl pads 2 Ns, so residue 0 is N N + 1 coding nt
        cds = {"T1": Seq("NN" + "A" * 10 + "C" * 6)}
        si = Superisoform("G", 1, {"T1": [(100, 109, "1"), (200, 205, "0")]}, cds)
        self.assertEqual(si.padding["T1"], 2)
        self.assertEqual(si.segments[1].phase, 0)
        self.assertEqual(si.transcript_nt_sequence("T1"), "A" * 10 + "C" * 6)
        self.assertEqual(si.residue_segments("T1"), [(0,), (0,), (0,), (0,), (1,), (1,)])

    def test_real_gene_rebuilds_cds(self):
        # Every HNF4A transcript's segments must concatenate back to its
        # Ensembl CDS sequence (minus padding and the separate stop codon).
        test_gene = gene.Gene.__new__(gene.Gene)
        test_gene.ensg_id = "ENSG00000101076"
        test_gene.strand = 1
        gtf_file = "reference_data/Homo_sapiens.GRCh38.109.gtf.gz"
        cds_file = "reference_data/Homo_sapiens.GRCh38.cds.all.fa.gz"
        si = test_gene.build_segment_superisoform(gtf_file, cds_file)
        cds = Transcript._get_cds_sequences(cds_file)
        self.assertGreater(len(si.transcript_segments), 0)
        for enst_id, blocks in si.transcript_blocks.items():
            pad = si.padding[enst_id]
            total = sum(end - start + 1 for start, end, _ in blocks)
            self.assertEqual(si.transcript_nt_sequence(enst_id), str(cds[enst_id][pad:pad + total]))


class TestDomainPositions(unittest.TestCase):
    def test_positions_are_one_based_inclusive(self):
        domain = Domain.from_positions_string("X", "104-352", "bc")
        self.assertEqual((int(domain.pos.start), int(domain.pos.end)), (103, 352))

    def test_lone_residue_is_one_residue_long(self):
        domain = Domain.from_positions_string("X", "171,246-247", "bc")
        self.assertEqual([len(part) for part in domain.pos.parts], [1, 2])


class TestSegmentMatching(unittest.TestCase):
    def test_split_codon_needs_consecutive_segments(self):
        # T1: 100-104 + 200-205; T2 inserts 150-152 (in frame) between them,
        # so the codon T1 splits across its two exons (bases 4-6) is joined
        # to a different exon in T2, while other shared codons still match.
        si = Superisoform("G", 1, {"T1": [(100, 104, "0"), (200, 205, "1")],
                                   "T2": [(100, 104, "0"), (150, 152, "1"), (200, 205, "1")]})
        codons = si.residue_codons("T1")
        self.assertEqual(len(codons), 3)
        self.assertTrue(si.contains_codon("T2", codons[0]))
        self.assertFalse(si.contains_codon("T2", codons[1]))
        self.assertTrue(si.contains_codon("T1", codons[1]))

    def test_choose_source_prefers_identity_then_isoform_1(self):
        xrefs = {"T1": [("P1", "Uniprot/SWISSPROT", 90.0)],
                 "T2": [("P1", "Uniprot/SWISSPROT", 100.0), ("P1-2", "Uniprot_isoform", None)],
                 "T3": [("P1", "Uniprot/SWISSPROT", 100.0), ("P1-1", "Uniprot_isoform", None)]}
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            sources = matching.choose_source_transcripts("G", ["T1", "T2", "T3"], xrefs)
        self.assertEqual(sources, {"P1": ("T3", 100.0)})

    def test_choose_source_falls_back_with_warning(self):
        xrefs = {"T1": [("P1", "Uniprot/SPTREMBL", 80.0)], "T2": [("P1", "Uniprot/SPTREMBL", 95.0)]}
        with self.assertWarns(UserWarning):
            sources = matching.choose_source_transcripts("G", ["T1", "T2"], xrefs)
        self.assertEqual(sources, {"P1": ("T2", 95.0)})

    def test_dpm1_one_row_per_site_per_transcript(self):
        # DPM1 (O60762) has 6 PPI binding sites; every protein-coding
        # transcript gets one row per site, and the source transcript
        # covers its own sites fully.
        reference_data = "reference_data/"
        test_gene = gene.Gene("ENSG00000000419", reference_data + "ppi_binding_sites.tsv",
                              reference_data + "Homo_sapiens.GRCh38.cds.all.fa.gz",
                              reference_data + "Homo_sapiens.GRCh38.109.uniprot.tsv.gz",
                              reference_data + "Homo_sapiens.GRCh38.109.gtf.gz",
                              reference_data + "Homo_sapiens.GRCh38.interpro_domains.tsv.gz",
                              domain_filter=["ppi_bs"], merge_overlapping_domains=False)
        rows = test_gene.coverage_rows
        self.assertEqual(len({row["domain_id"] for row in rows}), 6)
        self.assertEqual(len(rows), 6 * len(test_gene.superisoform.transcript_segments))
        own = [row for row in rows if row["transcript_id"] == row["source_transcript_id"]]
        self.assertEqual({row["source_transcript_id"] for row in own}, {"ENST00000371588"})
        self.assertTrue(all(row["coverage"] == 1.0 for row in own))


class TestGene(unittest.TestCase):
    def test_Gene_creation(self):
        test_ensg_id = "ENSG00000101076"
        ppi_binding_site_file = "/mnt/data/storage/WPI/Korkin_Lab/DNA_Binding_IF/reference_data/ppi_binding_sites.tsv"
        cds_fasta_file = "/mnt/data/storage/WPI/Korkin_Lab/DNA_Binding_IF/reference_data/Homo_sapiens.GRCh38.cds.all.fa.gz"
        uniprot_mapping_file = "/mnt/data/storage/WPI/Korkin_Lab/DNA_Binding_IF/reference_data/Homo_sapiens.GRCh38.109.uniprot.tsv.gz"
        gtf_file = "/mnt/data/storage/WPI/Korkin_Lab/DNA_Binding_IF/reference_data/Homo_sapiens.GRCh38.109.gtf.gz"
        interpro_domains_file = "/mnt/data/storage/WPI/Korkin_Lab/DNA_Binding_IF/reference_data/Homo_sapiens.GRCh38.interpro_domains.tsv.gz"
        test_gene = gene.Gene(test_ensg_id, ppi_binding_site_file, cds_fasta_file, uniprot_mapping_file,
                              gtf_file, interpro_domains_file, matchmode="alignment")
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
