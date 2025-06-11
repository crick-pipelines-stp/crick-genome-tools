"""
Tests covering the barcode_demux module
"""

# pylint: disable=missing-function-docstring,missing-class-docstring,no-member

import gzip
import os
import re

import pytest
from assertpy import assert_that

from crick_genome_tools.seq.barcode_demux import (
    assert_min_hamming_above_threshold,
    crosscheck_barcode_proximity,
    custom_priority_by_length_sort_key,
    demultiplex_fastq_by_barcode,
    extract_index_from_header_illumina,
    find_closest_match,
    find_min_hamming_distances,
    group_samples_by_index_length,
    hamming_distance,
    trim_merge_string,
    index_to_match_key,
    build_bk_tree_index
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
            "sample_5": {"key": "value", "index": "AAATA", "index2": "ATCTT"},
        }
        expected_dict = {
            (4, 0): {"sample_1": "ACGT"},
            (4, 4): {"sample_2": "ACGT,AGGT"},
            (8, 0): {"sample_3": "ACGTAGCT"},
        }
        expected_dict_2 = {
            (4, 0): {"sample_1": "ACGT"},
            (4, 4): {"sample_2": "ACGT+AGGT"},
            (8, 0): {"sample_3": "ACGTAGCT"},
            (5, 5): {"sample_4": "ACGTA:AGCTT", "sample_5": "AAATA,ATCTT"},
        }

        # Test and assert
        assert_that(group_samples_by_index_length(input_dict)).is_equal_to(expected_dict)
        assert_that(group_samples_by_index_length(input_dict_2)).is_equal_to(expected_dict_2)

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

    def test_custom_priority_by_length_sort_key_invalid(self):
        # Test and assert
        assert_that(custom_priority_by_length_sort_key).raises(TypeError).when_called_with(None)
        assert_that(custom_priority_by_length_sort_key).raises(TypeError).when_called_with("invalid_input")
        assert_that(custom_priority_by_length_sort_key).raises(TypeError).when_called_with([1234])

    def test_custom_priority_by_length_sort_key_isvalid_tuple_input(self):
        # Normal 2-tuple, no zero
        assert_that(custom_priority_by_length_sort_key((5, 5))).is_equal_to((False, -10))  # no 0 values, so returns False
        assert_that(custom_priority_by_length_sort_key((4, 0))).is_equal_to((True, -4))  # 0 value, so returns True
        assert_that(custom_priority_by_length_sort_key(6)).is_equal_to((True, -6))

        # Tie-breaker: same total length, but one has zero
        assert custom_priority_by_length_sort_key((4, 4)) < custom_priority_by_length_sort_key((8, 0))

    def test_custom_priority_by_length_sort_key_isvalid_dict_input(self):
        # Setup
        sample_barcode_dict = {
            (4, 0): {"sample_1": "ACGT"},
            (4, 4): {"sample_2": "ACGT+AGGT"},
            (8, 0): {"sample_3": "ACGTAGCT"},
            (5, 5): {"sample_4": "ACGTA:AGCTT"},
        }
        expected_sorted_samples = [(5, 5), (4, 4), (8, 0), (4, 0)]

        # Test and assert
        sorted_samples = sorted(sample_barcode_dict.keys(), key=custom_priority_by_length_sort_key)
        print(sorted_samples)

        assert_that(sorted_samples).is_equal_to(expected_sorted_samples)

    def test_trim_merge_string_isnone(self):
        assert_that(trim_merge_string).raises(ValueError).when_called_with(None, 3)
        assert_that(trim_merge_string).raises(ValueError).when_called_with("string", None)

    def test_trim_merge_string_invalid(self):
        # Test and assert
        assert_that(trim_merge_string).raises(ValueError).when_called_with("ACGTAGGT", -3)
        assert_that(trim_merge_string).raises(ValueError).when_called_with(["not_a_string"], 3)

    def test_trim_merge_string_isvalid(self):
        # Test and assert
        assert_that(trim_merge_string("ACGTAGGT", 5)).is_equal_to("ACGTA")
        assert_that(trim_merge_string("ACGTAGGT", 0)).is_equal_to("")
        assert_that(trim_merge_string("ACGT AGGT", 4)).is_equal_to("ACAG")
        assert_that(trim_merge_string("ACGT AGGT", 3)).is_equal_to("AA")
        assert_that(trim_merge_string("ACGT AGGT", 8)).is_equal_to("ACGTAGGT")
        assert_that(trim_merge_string("ACGTAAAC AGGT", 8)).is_equal_to("ACGTAGGT")
        assert_that(trim_merge_string("ACGT AGGT", 10)).is_equal_to("ACGTAGGT")

    def test_index_to_match_key_isvalid(self):
        # Setup
        name = "@LH00442:107:22YHM5LT3:2:1191:41040:4635 1:N:0:TTCAGAGGTC+CTGGAGCATC"
        sorted_grouped_samples_by_length = [
            (10, 10),
            (8, 0),
            (4, 0),
        ]
        grouped_samples_by_length = {
            (4, 0): {"sample_2": "ACGT"},
            (10, 10): {"sample_3": "TTCAGAGGTC,CTGGAGCATC"},
            (8, 0): {"sample_1": "ATTTAGAT"},
        }
        max_hamming_distance = 1

        for length, sample_dict in grouped_samples_by_length.items():
            print(sample_dict)
            grouped_bk_trees = build_bk_tree_index(sample_dict)

            # Test and assert
            match_key = index_to_match_key(name, sorted_grouped_samples_by_length, grouped_bk_trees, max_hamming_distance)
        assert_that(match_key).is_equal_to(("sample_3", "TTCAGAGGTC+CTGGAGCATC"))

    # def test_demultiplex_fastq_by_barcode_invalid_hamming_input(self, tmp_path):
    #     # Setup
    #     fastq_file = "tests/data/seq/L002_R1.fastq"
    #     barcode_sample_dict = {
    #         "sample_1": "ACTGGTGTCG-CAAGTCCTGT",
    #         "sample_2": "AGGTGGCTAC+CCACGTAACG",
    #         "sample_3": "TATCACTCTC+ACCTTGTTCT",
    #         "sample_4": "AGGTGGCTAC+CCACGTAACG",
    #     }
    #     max_hamming_distance = 1
    #     output_dir = tmp_path

    #     # Test and assert
    #     assert_that(demultiplex_fastq_by_barcode).raises(ValueError).when_called_with(
    #         barcode_sample_dict, fastq_file, max_hamming_distance, output_dir
    #     )
    #     demultiplex_fastq_by_barcode(fastq_file, barcode_sample_dict, max_hamming_distance, output_dir)

    @pytest.mark.parametrize(
        "fastq_file, barcode_sample_dict, max_hamming_distance, expected_samples, expected_file_content",
        [
            (
                "tests/data/seq/undetermined_L002_R1.fastq",
                {
                    "sample_1": "ACTT,NTAT",
                    "sample_2": "ACGT,AGGT",
                    "sample_3": "ACGTA",
                },
                1,
                ["sample_1", "sample_2", "sample_3", "undetermined"],
                {
                    "sample_1": "@LH00442:107:22YHM5LT3:2:1101:1000:1064 1:N:0:ACTTGACTAG+NTATCAACGG\nGNAGGGGCGGCCCGGCCCCCACCCCCACGCCCGCCCGGGAGGCGGACGGGGGGAGAGGGAGAGCGCGGCGACGGGTATCTGGCTTCCTCGGCCCCGGGATTCGGCGAAAGCTGCGGCCGGAGGGCTGTAACACTCGGGGTGAGGTGGTAGA\n+\nI#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIII9IIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",  # pytest: disable=line-too-long
                    "sample_2": "",
                    "sample_3": "",
                    "undetermined": "@LH00442:107:22YHM5LT3:2:1101:1092:1064 1:N:0:TCACCAGGAC+NCCTTGTCTC\nCNGCCACCTCCTCGGTCGCGCTGGCCGGGCCACCCGGGGTCAAAGCCACCTCACCCGAGCAAGTGGGTGCTAGTGAGGGCCGGGGGCGCCAGGCAGCACGGCAAGCGGAAGAGCCGAGCCGCAGCTCCGCAGCTGCCGGCGCCCGGGGAGA\n+\nI#IIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII*I9IIIIIIIII9IIIII9IIIII99I9III9II9IIIII999III99I99*9*I**9*II*99*\n@LH00442:107:22YHM5LT3:2:1101:1111:1064 1:N:0:GCGCTTCTAC+NTCCTTGGCT\nANTGACCTGTCATTTCAGCATGTCACCCCCAAGCCATCTCTAGGTGTACTTCTTCCATCGAGGAGAAAAATGTCTCTTTGACTTCTTAATGACACCGTGACGTTTGGTTCCAAAAAGGTGCCCTGGTAAATCTCCAGAAACACATTAGTTA\n+\nI#IIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIII*IIIIIIIIII9IIIII9IIIII99I9III9II9IIIII999III99I99*9*I**9*II*99*\n",  # pytest: disable=line-too-long
                },
            ),
            # (
            #     "tests/data/seq/undetermined_L002_R2.fastq",
            #     {"sample_A": "ACTTGACTAG+NTATCAACGG", "sample_B": "GCGCTTCTAC,NTCCTTGGCT", "sample_C": "ACCTTA+ACCTTA"},
            #     1,
            #     ["sample_A", "sample_B", "sample_C", "undetermined"],
            #     {
            #         "sample_A": "@LH00442:107:22YHM5LT3:2:1101:1000:1064 2:N:0:ACTTGACTAG+NTATCAACGG\nACCACCTCACCCCGAGTGTTACAGCCCTCCGGCCGCAGCTTTCGCCGAATCCCGGGGCCGAGGAAGCCAGATACCCGTCGCCGCGCTCTCCCTCTCCCCCCGTCCGCCTCCCGGGCGGGCGTGGGGGTGGGGGCCGGGCCGCCCCTCCAGA\n+\nIIII*IIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIII\n",  # pytest: disable=line-too-long
            #         "sample_B": "@LH00442:107:22YHM5LT3:2:1101:1111:1064 2:N:0:GCGCTTCTAC+NTCCTTGGCT\nAGACCTGCTGGGCTGACCACAGGCCTACAAACACGGACACTGCCTGAGAATAACTAATGTGTTTCTGGAGATTTACCAGGGCACCTTTTTGGAACCAAACGTCACGGTGTCATTACGAATTCAAAGAGACATCTTTCTCCTCGATGGAAGA\n+\nIIIIIIIIIIIII\n",
            #         "sample_C": "",
            #         "undetermined": "@LH00442:107:22YHM5LT3:2:1101:1092:1064 2:N:0:TCACCAGGAC+NCCTTGTCTC\nTGGAGACTCGCTGCCCGGGCGCCGGCAGCTGCGGAGCTGCGGCTCGGCTCTTCCGCTTGCCGTGCTGCCTGGCGCCCCCGGCCCTCACTAGCACCCACTTGCTCGGGTGAGGTGGCTTTGACCCCGGGTGGCCCGGCCAGCGCGACCGAGG\n+\nIIIIIIIIIIIIIIIIIII*IIIIIIIII9I9IIII9III9IIIIIIIII9III9III9IIII9999II99I9IIIIIIIIIII9IIIIIIIIIIIII9II9I99I9IIIIIIIIIIII9IIIII9IIIIIIII99II9IIIIIIIIII*9\n",  # pytest: disable=line-too-long
            #     },
            # ),
        ],
    )
    def test_demultiplex_fastq_by_barcode_valid(
        self, tmp_path, fastq_file, barcode_sample_dict, max_hamming_distance, expected_samples, expected_file_content
    ):  # pylint: disable=too-many-arguments
        # Setup
        output_dir = tmp_path
        # output_dir = "tests/data/seq/output"

        # Test
        demultiplex_fastq_by_barcode(barcode_sample_dict, fastq_file, max_hamming_distance, output_dir)

        # Assert
        read_lane_info = fastq_file.split("_", 1)[1].split(".", 1)[0]
        for sample in expected_samples:
            sample_file_path = os.path.join(output_dir, f"{sample}_{read_lane_info}.fastq.gz")

            # Check that the expected files were created
            assert_that(sample_file_path).exists()

            # Check that the content of each file is as expected
            with gzip.open(sample_file_path, "r") as f:
                file_content = f.read().decode("UTF-8")
                print(f"File content for {sample}: {file_content[:100]}...")  # Print first 100 characters for debugging

            assert_that(file_content).is_equal_to(expected_file_content[sample])
        # raise ValueError

    @pytest.mark.parametrize(
        "fastq_file, barcode_sample_dict, max_hamming_distance, expected_samples, expected_read_count",
        [
            (
                "tests/data/seq/sub_read_S1_L002_R1.fastq",
                {
                    "sample_1": "ACTT,NTAT",
                    "sample_2": "ACGT,AGGT",
                    "sample_3": "ACGTA",
                },
                1,
                ["sample_1", "sample_2", "sample_3", "undetermined"],
                {
                    "sample_1": 5,
                    # "sample_2": 0,
                    "sample_3": 57,
                    "undetermined": 938,
                },
            ),
        ],
    )
    def test_demultiplex_fastq_by_barcode_read_count_valid(
        self, tmp_path, fastq_file, barcode_sample_dict, max_hamming_distance, expected_samples, expected_read_count
    ):  # pytest: disable=line-too-long
        # Setup
        output_dir = tmp_path
        # output_dir = "tests/data/seq/output"

        # Test
        read_count = demultiplex_fastq_by_barcode(barcode_sample_dict, fastq_file, max_hamming_distance, output_dir)

        # Assert
        assert_that(read_count).is_equal_to(expected_read_count)

        read_lane_info = re.search(r"S\d+_([^\.]+)", fastq_file)
        if not read_lane_info:
            read_lane_info = fastq_file.split("_", 1)[1].split(".", 1)[0]
        for sample in expected_samples:
            sample_file_path = os.path.join(output_dir, f"{sample}_{read_lane_info}.fastq.gz")

            # Check that the expected files were created
            assert_that(sample_file_path).exists()

            # Check that the content of each file is as expected
            with gzip.open(sample_file_path, "r") as f:
                file_content = f.read().decode("UTF-8")
                newline_count = file_content.count("@")

            if sample in expected_read_count:
                assert_that(newline_count).is_equal_to(expected_read_count[sample])
            else:
                assert_that(newline_count).is_equal_to(0)

    # @pytest.mark.parametrize(
    #     "fastq_file_r1, fastq_file_r2, barcode_sample_dict, max_hamming_distance, expected_samples, expected_file_content_r1, expected_file_content_r2, expected_read_count",
    #     [
    #         (
    #             "tests/data/seq/undetermined_L002_R1.fastq",
    #             "tests/data/seq/undetermined_L002_R2.fastq",
    #             {
    #                 "sample_1": "ACTT,NTAT",
    #                 "sample_2": "ACGT,AGGT",
    #                 "sample_3": "ACGTA",
    #             },
    #             1,
    #             ["sample_1", "sample_2", "sample_3", "undetermined"],
    #             {
    #                 "sample_1": "@LH00442:107:22YHM5LT3:2:1101:1000:1064 1:N:0:ACTTGACTAG+NTATCAACGG\nGNAGGGGCGGCCCGGCCCCCACCCCCACGCCCGCCCGGGAGGCGGACGGGGGGAGAGGGAGAGCGCGGCGACGGGTATCTGGCTTCCTCGGCCCCGGGATTCGGCGAAAGCTGCGGCCGGAGGGCTGTAACACTCGGGGTGAGGTGGTAGA\n+\nI#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIII9IIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",  # pytest: disable=line-too-long
    #                 "sample_2": "",
    #                 "sample_3": "",
    #                 "undetermined": "@LH00442:107:22YHM5LT3:2:1101:1092:1064 1:N:0:TCACCAGGAC+NCCTTGTCTC\nCNGCCACCTCCTCGGTCGCGCTGGCCGGGCCACCCGGGGTCAAAGCCACCTCACCCGAGCAAGTGGGTGCTAGTGAGGGCCGGGGGCGCCAGGCAGCACGGCAAGCGGAAGAGCCGAGCCGCAGCTCCGCAGCTGCCGGCGCCCGGGGAGA\n+\nI#IIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII*I9IIIIIIIII9IIIII9IIIII99I9III9II9IIIII999III99I99*9*I**9*II*99*\n@LH00442:107:22YHM5LT3:2:1101:1111:1064 1:N:0:GCGCTTCTAC+NTCCTTGGCT\nANTGACCTGTCATTTCAGCATGTCACCCCCAAGCCATCTCTAGGTGTACTTCTTCCATCGAGGAGAAAAATGTCTCTTTGACTTCTTAATGACACCGTGACGTTTGGTTCCAAAAAGGTGCCCTGGTAAATCTCCAGAAACACATTAGTTA\n+\nI#IIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIII*IIIIIIIIII9IIIII9IIIII99I9III9II9IIIII999III99I99*9*I**9*II*99*\n",  # pytest: disable=line-too-long
    #             },
    #             {
    #                 "sample_1": "@LH00442:107:22YHM5LT3:2:1101:1000:1064 2:N:0:ACTTGACTAG+NTATCAACGG\nACCACCTCACCCCGAGTGTTACAGCCCTCCGGCCGCAGCTTTCGCCGAATCCCGGGGCCGAGGAAGCCAGATACCCGTCGCCGCGCTCTCCCTCTCCCCCCGTCCGCCTCCCGGGCGGGCGTGGGGGTGGGGGCCGGGCCGCCCCTCCAGA\n+\nIIII*IIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIII\n",  # pytest: disable=line-too-long
    #                 "sample_2": "",
    #                 "sample_3": "",
    #                 "undetermined": "@LH00442:107:22YHM5LT3:2:1101:1092:1064 2:N:0:TCACCAGGAC+NCCTTGTCTC\nTGGAGACTCGCTGCCCGGGCGCCGGCAGCTGCGGAGCTGCGGCTCGGCTCTTCCGCTTGCCGTGCTGCCTGGCGCCCCCGGCCCTCACTAGCACCCACTTGCTCGGGTGAGGTGGCTTTGACCCCGGGTGGCCCGGCCAGCGCGACCGAGG\n+\nIIIIIIIIIIIIIIIIIII*IIIIIIIII9I9IIII9III9IIIIIIIII9III9III9IIII9999II99I9IIIIIIIIIII9IIIIIIIIIIIII9II9I99I9IIIIIIIIIIII9IIIII9IIIIIIII99II9IIIIIIIIII*9\n@LH00442:107:22YHM5LT3:2:1101:1111:1064 2:N:0:GCGCTTCTAC+NTCCTTGGCT\nAGACCTGCTGGGCTGACCACAGGCCTACAAACACGGACACTGCCTGAGAATAACTAATGTGTTTCTGGAGATTTACCAGGGCACCTTTTTGGAACCAAACGTCACGGTGTCATTACGAATTCAAAGAGACATCTTTCTCCTCGATGGAAGA\n+\nIIIIIIIIIIIII\n",  # pytest: disable=line-too-long
    #             },
    #             {
    #                 "sample_1": 1,
    #                 # "sample_2": 0,
    #                 # "sample_3": 0,
    #                 "undetermined": 2,
    #             },
    #         ),
    #     ],
    # )
    # def test_demultiplex_fastq_by_barcode_valid_dualindex(
    #     self,
    #     tmp_path,
    #     fastq_file_r1,
    #     fastq_file_r2,
    #     barcode_sample_dict,
    #     max_hamming_distance,
    #     expected_samples,
    #     expected_file_content_r1,
    #     expected_file_content_r2,
    #     expected_read_count,
    # ):  # pylint: disable=too-many-arguments
    #     # Setup
    #     output_dir = tmp_path
    #     # output_dir = "tests/data/seq/output"

    #     # Test
    #     demultiplex_fastq_by_barcode(barcode_sample_dict, fastq_file_r1, max_hamming_distance, output_dir, fastq_file_r2)

    #     # Assert
    #     for sample in expected_samples:
    #         sample_file_path_r1 = os.path.join(output_dir, f"{sample}_L002_R1.fastq.gz")
    #         sample_file_path_r2 = os.path.join(output_dir, f"{sample}_L002_R2.fastq.gz")

    #         # Check that the expected files were created
    #         assert_that(sample_file_path_r1).exists()
    #         assert_that(sample_file_path_r2).exists()

    #         # Check that the content of each file is as expected
    #         with gzip.open(sample_file_path_r1, "r") as f:
    #             file_content_r1 = f.read().decode("UTF-8")
    #             newline_count_r1 = file_content_r1.count("@")

    #         with gzip.open(sample_file_path_r2, "r") as f:
    #             file_content_r2 = f.read().decode("UTF-8")
    #             newline_count_r2 = file_content_r2.count("@")

    #         assert_that(file_content_r1).is_equal_to(expected_file_content_r1[sample])
    #         assert_that(file_content_r2).is_equal_to(expected_file_content_r2[sample])

    #         # Check the read count in each file is as expected
    #         if sample in expected_read_count:
    #             assert_that(newline_count_r1).is_equal_to(expected_read_count[sample])
    #         else:
    #             assert_that(newline_count_r1).is_equal_to(0)

    #         if sample in expected_read_count:
    #             assert_that(newline_count_r2).is_equal_to(expected_read_count[sample])
    #         else:
    #             assert_that(newline_count_r2).is_equal_to(0)
