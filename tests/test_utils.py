from unittest import TestCase
from scripts.utils import reverse_complement

class TestUtils(TestCase):
    def test_reverse_complement_1(self):
        seq = "A"
        actual = reverse_complement(seq)
        self.assertEqual(actual, "T")

    def test_reverse_complement_2(self):
        seq = "GGCGCTAAAAATAGCGACTTGGGCGATTTTTGCAGCAAACGATTCAAAAGATGA"
        actual = reverse_complement(seq)
        self.assertEqual(actual, "TCATCTTTTGAATCGTTTGCTGCAAAAATCGCCCAAGTCGCTATTTTTAGCGCC")

    def test_reverse_complement_3(self):
        seq = "CGCCTTTGTTCCAGCCGCCTTCGACCTGATGAGCGGCAACTGCGCCGCCCCATAAGAAATCTT"
        actual = reverse_complement(seq)
        self.assertEqual(actual, "AAGATTTCTTATGGGGCGGCGCAGTTGCCGCTCATCAGGTCGAAGGCGGCTGGAACAAAGGCG")

