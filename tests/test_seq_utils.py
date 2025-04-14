"""
Tests for seq utils.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring


import pytest
from assertpy import assert_that

from crick_genome_tools.seq.utils import cumulative_hamming_distance, hamming_distance, rev_comp


class TestHammingDistance:

    def test_seq_util_hamming_distance_same_length(self):
        assert_that(hamming_distance("ACTG", "ACTG")).is_equal_to(0)
        assert_that(hamming_distance("ACTG", "ACCG")).is_equal_to(1)
        assert_that(hamming_distance("ACTG", "ACCA")).is_equal_to(2)
        assert_that(hamming_distance("ACTG", "TGCA")).is_equal_to(4)

    def test_seq_util_hamming_distance_different_length(self):
        assert_that(hamming_distance("ACTG", "ACT")).is_equal_to(1)
        assert_that(hamming_distance("ACT", "ACTG")).is_equal_to(1)
        assert_that(hamming_distance("ACTG", "AC")).is_equal_to(2)
        assert_that(hamming_distance("AC", "ACTG")).is_equal_to(2)

    def test_seq_util_hamming_distance_empty_sequences(self):
        assert_that(hamming_distance("", "")).is_equal_to(0)
        assert_that(hamming_distance("ACTG", "")).is_equal_to(4)
        assert_that(hamming_distance("", "ACTG")).is_equal_to(4)

    def test_seq_util_hamming_distance_mixed_cases(self):
        assert_that(hamming_distance("ACTG", "actg")).is_equal_to(4)
        assert_that(hamming_distance("aCtG", "AcTg")).is_equal_to(4)

    def test_seq_util_hamming_distance_special_characters(self):
        assert_that(hamming_distance("A!C@G", "A!C@G")).is_equal_to(0)
        assert_that(hamming_distance("A!C@G", "A!C@T")).is_equal_to(1)
        assert_that(hamming_distance("A!C@G", "A!C@")).is_equal_to(1)
        assert_that(hamming_distance("A!C@", "A!C@G")).is_equal_to(1)

    def test_seq_util_hamming_distance_numeric_characters(self):
        assert_that(hamming_distance("1234", "1234")).is_equal_to(0)
        assert_that(hamming_distance("1234", "1235")).is_equal_to(1)
        assert_that(hamming_distance("1234", "124")).is_equal_to(2)
        assert_that(hamming_distance("124", "1234")).is_equal_to(2)

    def test_seq_util_hamming_distance_unicode_characters(self):
        assert_that(hamming_distance("你好", "你好")).is_equal_to(0)
        assert_that(hamming_distance("你好", "你坏")).is_equal_to(1)
        assert_that(hamming_distance("你好", "你")).is_equal_to(1)
        assert_that(hamming_distance("你", "你好")).is_equal_to(1)

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
    def test_seq_util_cumulative_hamming_distance(self, dict1, dict2, expected):
        assert cumulative_hamming_distance(dict1, dict2) == expected

    def test_seq_util_rev_comp_basic(self):
        assert_that(rev_comp("ATGC")).is_equal_to("GCAT")
        assert_that(rev_comp("AATT")).is_equal_to("AATT")
        assert_that(rev_comp("CCGG")).is_equal_to("CCGG")

    def test_seq_util_rev_comp_empty(self):
        assert_that(rev_comp("")).is_equal_to("")

    def test_seq_util_rev_comp_single_nucleotide(self):
        assert_that(rev_comp("A")).is_equal_to("T")
        assert_that(rev_comp("C")).is_equal_to("G")

    def test_seq_util_rev_comp_mixed_case(self):
        assert_that(rev_comp).raises(KeyError).when_called_with("aTgc")

    def test_seq_util_rev_comp_invalid_characters(self):
        assert_that(rev_comp("ATGCNNN")).is_equal_to("NNNGCAT")
