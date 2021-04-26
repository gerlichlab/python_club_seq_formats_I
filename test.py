import unittest
from Bio import SeqIO
import lib.homework as hw


class TestGetNucleotideContent(unittest.TestCase):

    def test_example_w_one_sequence(self):
        """Tests example fasta file with one sequence"""
        with open("./example_files/example2.fa") as f:
            fasta = SeqIO.parse(f, "fasta")
            # run function
            result = hw.get_nucleotide_content(fasta)
        # compare result
        expected = {'HSGLTH1': {'A': 0.14313725490196078,
                                'C': 0.35490196078431374,
                                'G': 0.35,
                                'T': 0.15196078431372548}}
        # test names equal
        self.assertEqual(expected.keys(), result.keys())
        # test values equal -> are floats, need to test almostequal
        for key in expected:
            for nuc in expected[key]:
                self.assertAlmostEqual(result[key][nuc], expected[key][nuc])

    def test_example_w_two_sequences(self):
        """Tests example fasta file with one sequence"""
        with open("./example_files/example1.fa") as f:
            fasta = SeqIO.parse(f, "fasta")
            # run function
            result = hw.get_nucleotide_content(fasta)
        # compare result
        expected = {'HSBGPG': {'A': 0.17627944760357434,
                               'C': 0.3322502030869212,
                               'G': 0.30300568643379366,
                               'T': 0.18846466287571081},
                    'HSGLTH1': {'A': 0.14313725490196078,
                                'C': 0.35490196078431374,
                                'G': 0.35,
                                'T': 0.15196078431372548}}
        # test names equal
        self.assertEqual(expected.keys(), result.keys())
        # test values equal -> are floats, need to test almostequal
        for key in expected:
            for nuc in expected[key]:
                self.assertAlmostEqual(result[key][nuc], expected[key][nuc])


if __name__ == "__main__":
    res = unittest.main(verbosity=3, exit=False)