"""
Tests for seq utils.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring


import unittest

from crick_genome_tools.seq.utils import hamming_distance

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
