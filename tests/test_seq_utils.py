"""
Tests for seq utils.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring


import unittest

import pytest

from crick_genome_tools.seq.utils import cumulative_hamming_distance, hamming_distance


class TestHammingDistance(unittest.TestCase):

    def test_hamming_distance_same_length(self):
        self.assertEqual(hamming_distance("ACTG", "ACTG"), 0)
        self.assertEqual(hamming_distance("ACTG", "ACCG"), 1)
        self.assertEqual(hamming_distance("ACTG", "ACCA"), 2)
        self.assertEqual(hamming_distance("ACTG", "TGCA"), 4)

    def test_hamming_distance_different_length(self):
        self.assertEqual(hamming_distance("ACTG", "ACT"), 1)
        self.assertEqual(hamming_distance("ACT", "ACTG"), 1)
        self.assertEqual(hamming_distance("ACTG", "AC"), 2)
        self.assertEqual(hamming_distance("AC", "ACTG"), 2)

    def test_hamming_distance_empty_sequences(self):
        self.assertEqual(hamming_distance("", ""), 0)
        self.assertEqual(hamming_distance("ACTG", ""), 4)
        self.assertEqual(hamming_distance("", "ACTG"), 4)

    def test_hamming_distance_mixed_cases(self):
        self.assertEqual(hamming_distance("ACTG", "actg"), 4)
        self.assertEqual(hamming_distance("aCtG", "AcTg"), 4)

    def test_hamming_distance_special_characters(self):
        self.assertEqual(hamming_distance("A!C@G", "A!C@G"), 0)
        self.assertEqual(hamming_distance("A!C@G", "A!C@T"), 1)
        self.assertEqual(hamming_distance("A!C@G", "A!C@"), 1)
        self.assertEqual(hamming_distance("A!C@", "A!C@G"), 1)

    def test_hamming_distance_numeric_characters(self):
        self.assertEqual(hamming_distance("1234", "1234"), 0)
        self.assertEqual(hamming_distance("1234", "1235"), 1)
        self.assertEqual(hamming_distance("1234", "124"), 2)
        self.assertEqual(hamming_distance("124", "1234"), 2)

    def test_hamming_distance_unicode_characters(self):
        self.assertEqual(hamming_distance("你好", "你好"), 0)
        self.assertEqual(hamming_distance("你好", "你坏"), 1)
        self.assertEqual(hamming_distance("你好", "你"), 1)
        self.assertEqual(hamming_distance("你", "你好"), 1)


class TestHammingDistanceWithFixtures:

    @pytest.mark.parametrize(
        "dict1, dict2, expected",
        [
            # Both dictionaries are empty, so the cumulative Hamming distance is 0
            ({}, {}, 0),
            # Both dictionaries are identical, so the cumulative Hamming distance is 0
            ({"seq1": "AGCT", "seq2": "TCGA"}, {"seq1": "AGCT", "seq2": "TCGA"}, 0),
            # Both dictionaries have sequences with 1 mismatch each, total distance is 2
            ({"seq1": "AGCT", "seq2": "TCGA"}, {"seq1": "AGTT", "seq2": "TCCA"}, 2),
            # Each dictionary has an extra sequence, total distance is 8
            ({"seq1": "AGCT", "seq2": "TCGA"}, {"seq1": "AGCT", "seq3": "TGCA"}, 8),
            # Second dictionary is empty, total distance is the sum of lengths of sequences in the first dictionary
            ({"seq1": "AGCT", "seq2": "TCGA"}, {}, 8),
            # Sequences have different lengths, distance is the difference in lengths plus mismatches
            ({"seq1": "AGCT", "seq2": "TCGA"}, {"seq1": "AGC", "seq2": "TCGAT"}, 2),
            # Single sequence with 2 mismatches
            ({"seq1": "AGCT"}, {"seq1": "AGTG"}, 2),
            # Sequences have different lengths and mismatches, total distance is 3
            ({"seq1": "TGCT", "seq2": "TCGA"}, {"seq1": "AGCTT", "seq2": "TCG"}, 3),
        ],
    )
    def test_cumulative_hamming_distance(self, dict1, dict2, expected):
        assert cumulative_hamming_distance(dict1, dict2) == expected
