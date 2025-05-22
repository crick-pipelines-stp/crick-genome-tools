"""
Tests covering the barcode_demux module
"""

# pylint: disable=missing-function-docstring,missing-class-docstring,no-member

import os

import pytest
from assertpy import assert_that

from crick_genome_tools.seq.barcode_demux import (
    assert_min_hamming_above_threshold,
    crosscheck_barcode_proximity,
    demultiplex_fastq_by_barcode,
    extract_index_from_header_illumina,
    find_closest_match,
    find_min_hamming_distances,
    find_sample_for_read_index,
    group_samples_by_index_length,
    hamming_distance,
    trim_merge_string,
)


class TestBarcodeDemux:
    def test_extract_index_from_header_none(self):
        assert_that(extract_index_from_header_illumina).raises(ValueError).when_called_with(None)

    def test_extract_index_from_header_empty(self):
        assert_that(extract_index_from_header_illumina).raises(ValueError).when_called_with("")

    def test_extract_index_from_header_invalid(self):
        header = "@LH00442:107:22YHM5LT3:2:1101:1092:1064 2:N:0:TCACCA:GGAC+NCCTTGT:CTC"
        assert_that(extract_index_from_header_illumina).raises(ValueError).when_called_with(header)

    def test_extract_index_from_header_isvalid(self):
        header = "@LH00442:107:22YHM5LT3:2:1101:1092:1064 2:N:0:TCACCAGGAC+NCCTTGTCTC"
        header_2 = "@LH00442:107 1:N:0:GCGCTTCTAC+NTCCTTGGCT"
        header_3 = "invalid_header"

        # Test and assert
        assert_that(extract_index_from_header_illumina(header)).is_equal_to("TCACCAGGAC+NCCTTGTCTC")
        assert_that(extract_index_from_header_illumina(header_2)).is_equal_to("GCGCTTCTAC+NTCCTTGGCT")
        assert_that(extract_index_from_header_illumina(header_3)).is_equal_to("")

    def test_group_samples_by_index_length_none(self):
        assert_that(group_samples_by_index_length).raises(TypeError).when_called_with(None)

    def test_group_samples_by_index_length_emptydict(self):
        assert_that(group_samples_by_index_length({})).is_equal_to({})

    def test_group_samples_by_index_length_invalid(self):
        # Call the fuction with a dictionary with invalid barcode values
        assert_that(group_samples_by_index_length).raises(TypeError).when_called_with({"sample_1": 1234})

        # Call the function with an input that is not a dict
        assert_that(group_samples_by_index_length).raises(TypeError).when_called_with("not_a_dict")

    def test_group_samples_by_index_length_valid(self):
        input_dict = {
            "sample_1": "ACGT",
            "sample_2": "ACGT,AGGT",
            "sample_3": "ACGTAGCT",
        }
        input_dict_2 = {
            "sample_1": "ACGT",
            "sample_2": "ACGT+AGGT",
            "sample_3": {"key": "value", "barcode": "ACGTAGCT"},
            "sample_4": {"key": "value", "index": "ACGTA:AGCTT"},
        }
        expected_dict = {
            4: {"sample_1": "ACGT"},
            8: {"sample_2": "ACGT,AGGT", "sample_3": "ACGTAGCT"},
        }
        expected_dict_2 = {
            4: {"sample_1": "ACGT"},
            8: {"sample_2": "ACGT+AGGT", "sample_3": "ACGTAGCT"},
            10: {"sample_4": "ACGTA:AGCTT"},
        }

        # Test and assert
        assert_that(group_samples_by_index_length(input_dict)).is_equal_to(expected_dict)
        assert_that(group_samples_by_index_length(input_dict_2)).is_equal_to(expected_dict_2)

    def test_find_sample_for_read_index_none(self):
        assert_that(find_sample_for_read_index).raises(ValueError).when_called_with(None, {"dict": "test"})

        assert_that(find_sample_for_read_index).raises(ValueError).when_called_with("str_input", None)

    def test_find_sample_for_read_index_invalid(self):
        assert_that(find_sample_for_read_index).raises(ValueError).when_called_with(["not_a_string_input"], {"dict": "test"})

        assert_that(find_sample_for_read_index).raises(ValueError).when_called_with("str_input", ["not_a_dict"])

    def test_find_sample_for_read_index_valid(self):
        # Setup
        read_index = "ACGT"
        barcode_dict = {
            4: {
                "sample_1": {"ACGT", "AGGT"},
                "sample_2": {"CGTA", "TGCA"},
            }
        }

        # Test and assert
        assert_that(find_sample_for_read_index(read_index, barcode_dict)).is_equal_to("sample_1")

    @pytest.mark.parametrize(
        "sequence1, sequence2",
        [
            ("AATCG", None),
            (None, "AATCG"),
        ],
    )
    def test_hamming_distance_isnone(self, sequence1, sequence2):
        assert_that(hamming_distance).raises(ValueError).when_called_with(sequence1, sequence2)

    @pytest.mark.parametrize(
        "sequence1, sequence2, expected_result",
        [
            ("ACGT", "ACGT", 0),
            ("ACGT", "AGGT", 1),
            ("ACGT", "ACGTT", 0),
            ("ACGT", "ACG", 0),
            ("ACGT", "invalid", 4),
            ("invalid", "ACGT", 4),
            ("ACCT", "ATAC", 3),
        ],
    )
    def test_hamming_distance_isvalid(self, sequence1, sequence2, expected_result):
        # Test and assert
        result = hamming_distance(sequence1, sequence2)
        assert_that(result).is_equal_to(expected_result)

    def test_find_closest_match_invalid_input(self):
        # Test and assert
        # not a dictionary
        assert_that(find_closest_match).raises(ValueError).when_called_with("ACGT", "ACGTTT", 1)
        # not a string
        assert_that(find_closest_match).raises(ValueError).when_called_with({"sample_1": "ACGT"}, 1234, 1)
        # not a number
        assert_that(find_closest_match).raises(ValueError).when_called_with({"sample_1": "ACGT"}, "ACGTTT", "not_a_number")

    @pytest.mark.parametrize(
        "sample_barcode_dict, sequence, max_hamming_distance, expected_result",
        [
            (
                {
                    "sample_1": "ACGTAGGT",
                    "sample_2": "ACGTAAAA",
                    "sample_3": "AAGTAGGG",
                },
                "AAGTAGGG",
                1,
                "sample_3",
            ),
            (
                {
                    "sample_1": "ACGTAGGT",
                    "sample_2": "ACGTAAAA",
                    "sample_3": "AAGTAGGG",
                },
                "AAGTAGGT",
                1,
                "sample_1",
            ),
            ({"sample_1": "ACGTAGGT"}, "AAAAAAA", 3, "undetermined"),
            (
                {
                    "sample_1": "TTTTTTTT",
                    "sample_2": "ACGTAGCT",
                },
                "ACCGAGCA",
                4,
                "sample_2",
            ),
        ],
    )
    def test_find_closest_match_isvalid(self, sample_barcode_dict, sequence, max_hamming_distance, expected_result):
        # Test and assert
        result = find_closest_match(sample_barcode_dict, sequence, max_hamming_distance)
        assert_that(result).is_equal_to(expected_result)

    def test_crosscheck_barcode_proximity_invalid_input(self):
        # Setup
        sample_barcode_list = ["invalid", "input", "not", "a", "dict"]

        # Test and assert
        assert_that(crosscheck_barcode_proximity).raises(ValueError).when_called_with(sample_barcode_list)
        assert_that(crosscheck_barcode_proximity).raises(ValueError).when_called_with("invalid_input")
        assert_that(crosscheck_barcode_proximity).raises(ValueError).when_called_with(1234)
        assert_that(crosscheck_barcode_proximity).raises(ValueError).when_called_with(None)

    def test_crosscheck_barcode_proximity_invalid_barcode_length(self):
        # Setup
        sample_barcode_dict = {
            "sample_1": "ACGTAGGT",
            "sample_2": "ACGTAA",
            "sample_3": "ACGGA",
        }

        # Test and assert
        assert_that(crosscheck_barcode_proximity).raises(ValueError).when_called_with(sample_barcode_dict)

    def test_crosscheck_barcode_proximity_single_barcode(self):
        # Setup
        sample_barcode_dict = {"sample_1": "ACGTAGGT"}

        # Test and assert
        assert_that(crosscheck_barcode_proximity(sample_barcode_dict)).is_equal_to([])
        # No pairs to compare, so the result should be an empty list

    def test_crosscheck_barcode_proximity_isvalid(self):
        # Setup
        sample_barcode_dict = {
            "sample_1": "ACGTAGGT",
            "sample_2": "ACGTAAAA",
            "sample_3": "ACGTAGGA",
        }
        expected_result = [("ACGTAGGT", "ACGTAAAA", 3), ("ACGTAGGT", "ACGTAGGA", 1), ("ACGTAAAA", "ACGTAGGA", 2)]
        # ('sample_1', 'sample_2', 3), ('sample_1', 'sample_3', 1), ('sample_2', 'sample_3', 2)

        # Test and assert
        assert_that(crosscheck_barcode_proximity(sample_barcode_dict)).is_equal_to(expected_result)

    def test_find_min_hamming_distances_invalid_input(self):
        # Test and assert
        # not a dictionary
        assert_that(find_min_hamming_distances).raises(ValueError).when_called_with("ACGT")  # string
        assert_that(find_min_hamming_distances).raises(ValueError).when_called_with(["sample_1"])  # list
        assert_that(find_min_hamming_distances).raises(ValueError).when_called_with(1)  # int
        assert_that(find_min_hamming_distances).raises(ValueError).when_called_with(None)  # None

    def test_find_min_hamming_distances_isvalid(self):
        # Setup
        sample_barcode_dict = {
            8: [("ACGTAGGT", "ACGTAAAA", 3), ("ACGTAAAT", "ACGTAGGA", 3), ("ACGTAAAA", "ACGTAGGA", 2)],
            4: [("ACGT", "ACAA", 3), ("ACGT", "ACGA", 1), ("AAAA", "AGGA", 2)],
        }

        # Test and assert
        assert_that(find_min_hamming_distances(sample_barcode_dict)).is_equal_to({8: 2, 4: 1})

    def test_assert_min_hamming_above_threshold_hamming_below_max_threshold(self):
        # Setup
        sample_barcode_dict = {8: 4, 4: 5, 3: 3}

        # Test and assert
        assert_that(assert_min_hamming_above_threshold).raises(ValueError).when_called_with(sample_barcode_dict, 4)
        
    def test_assert_min_hamming_above_threshold_isinvalid(self):
        # Setup
        sample_barcode_dict = {8: 4, 4: 5, 3: 3}
        assert_that(assert_min_hamming_above_threshold).raises(ValueError).when_called_with(sample_barcode_dict, None)
        assert_that(assert_min_hamming_above_threshold).raises(ValueError).when_called_with("invalid_input", 3)
        assert_that(assert_min_hamming_above_threshold).raises(ValueError).when_called_with(None, 3)

    def test_assert_min_hamming_above_threshold_isvalid(self):
        # Setup
        sample_barcode_dict = {8: 4, 4: 5, 3: 3}

        # Test and assert
        try:
            assert_min_hamming_above_threshold(sample_barcode_dict, 3)  # function that doesn't return anything
        except ValueError as e:
            assert False, f"Function raised an unexpected exception: {e}"
            # assertion can't be done using assertpy

    def test_trim_merge_string_isnone(self):
        assert_that(trim_merge_string).raises(ValueError).when_called_with(None, 3)
        assert_that(trim_merge_string).raises(ValueError).when_called_with("string", None)

    def test_trim_merge_string_invalid(self):
        # Test and assert
        assert_that(trim_merge_string).raises(ValueError).when_called_with("ACGTAGGT", -3)
        assert_that(trim_merge_string).raises(ValueError).when_called_with(["not_a_string"], 3)

    def test_trim_merge_string_isvalid(self):
        # Test and assert
        assert_that(trim_merge_string("ACGTAGGT", 3)).is_equal_to("ACGTA")
        assert_that(trim_merge_string("ACGTAGGT", 0)).is_equal_to("ACGTAGGT")
        assert_that(trim_merge_string("ACGT AGGT", 2)).is_equal_to("ACGAGG")
        assert_that(trim_merge_string("ACGT AGGT", 3)).is_equal_to("ACGAGG")
        assert_that(trim_merge_string("ACGT AGGT", 0)).is_equal_to("ACGTAGGT")
        assert_that(trim_merge_string("ACGTAAAC AGGT", 4)).is_equal_to("ACGTAAAG")

    @pytest.mark.parametrize(
        "fastq_file, barcode_sample_dict, max_hamming_distance, expected_samples, expected_file_content",
        [
            (
                "tests/data/seq/L002_R1.fastq",
                {
                    "sample_1": "ACTT,NTAT",
                    "sample_2": "ACGT,AGGT",
                    "sample_3": "ACGTA",
                },
                1,
                ["sample_1", "sample_2", "undetermined"],
                {
                    "sample_1": "GNAGGGGCGGCCCGGCCCCCACCCCCACGCCCGCCCGGGAGGCGGACGGGGGGAGAGGGAGAGCGCGGCGACGGGTATCTGGCTTCCTCGGCCCCGGGATTCGGCGAAAGCTGCGGCCGGAGGGCTGTAACACTCGGGGTGAGGTGGTAGA",
                    "sample_2": "",
                    "sample_3": "",
                    "undetermined": "CNGCCACCTCCTCGGTCGCGCTGGCCGGGCCACCCGGGGTCAAAGCCACCTCACCCGAGCAAGTGGGTGCTAGTGAGGGCCGGGGGCGCCAGGCAGCACGGCAAGCGGAAGAGCCGAGCCGCAGCTCCGCAGCTGCCGGCGCCCGGGGAGA\nANTGACCTGTCATTTCAGCATGTCACCCCCAAGCCATCTCTAGGTGTACTTCTTCCATCGAGGAGAAAAATGTCTCTTTGACTTCTTAATGACACCGTGACGTTTGGTTCCAAAAAGGTGCCCTGGTAAATCTCCAGAAACACATTAGTTA",
                },
            ),
            (
                "tests/data/seq/L002_R2.fastq",
                {
                    "sample_A": "ACTTGACTAG+NTATCAACGG",
                    "sample_B": "GCGCTTCTAC,NTCCTTGGCT",
                },
                1,
                ["sample_A", "sample_B", "undetermined"],
                {
                    "sample_A": "ACCACCTCACCCCGAGTGTTACAGCCCTCCGGCCGCAGCTTTCGCCGAATCCCGGGGCCGAGGAAGCCAGATACCCGTCGCCGCGCTCTCCCTCTCCCCCCGTCCGCCTCCCGGGCGGGCGTGGGGGTGGGGGCCGGGCCGCCCCTCCAGA",
                    "sample_B": "AGACCTGCTGGGCTGACCACAGGCCTACAAACACGGACACTGCCTGAGAATAACTAATGTGTTTCTGGAGATTTACCAGGGCACCTTTTTGGAACCAAACGTCACGGTGTCATTACGAATTCAAAGAGACATCTTTCTCCTCGATGGAAGA",
                    "undetermined": "TGGAGACTCGCTGCCCGGGCGCCGGCAGCTGCGGAGCTGCGGCTCGGCTCTTCCGCTTGCCGTGCTGCCTGGCGCCCCCGGCCCTCACTAGCACCCACTTGCTCGGGTGAGGTGGCTTTGACCCCGGGTGGCCCGGCCAGCGCGACCGAGG",
                },
            ),
        ],
    )
    def test_demultiplex_fastq_by_barcode_valid(
        self, tmp_path, fastq_file, barcode_sample_dict, max_hamming_distance, expected_samples, expected_file_content
    ):
        # Setup
        output_dir = tmp_path
        # output_dir = "tests/data/seq/output"

        # Test
        demultiplex_fastq_by_barcode(fastq_file, barcode_sample_dict, max_hamming_distance, output_dir)

        # Assert

        for sample in expected_samples:
            sample_file_path = os.path.join(output_dir, f"{sample}.txt")

            # Check that the expected files were created
            assert_that(sample_file_path).exists()

            # Check that the content of each file is as expected
            with open(sample_file_path, "r", encoding="utf-8") as f:
                file_content = f.read().strip()

            assert_that(file_content).is_equal_to(expected_file_content[sample])
            raise ValueError
