"""
Tests covering the barcode_demux module
"""

# pylint: disable=missing-function-docstring,missing-class-docstring,no-member

import os

import pytest
from assertpy import assert_that

from crick_genome_tools.seq.barcode_demux import (
    demultiplex_fastq_by_barcode,
    extract_index_from_header_illumina,
    find_sample_for_read_index,
    gen_nearby_seqs,
    generate_nearby_barcodes_by_length,
    group_samples_by_index_length,
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

    def test_gen_nearby_seqs_with_nodist(self):
        # Setup
        seq = "atcg"
        barcode_sets = [seq]

        # Test and Assert
        assert_that(list(gen_nearby_seqs(seq, barcode_sets))).is_equal_to(["atcg"])

    @pytest.mark.parametrize("maxdist", [1, 2])
    @pytest.mark.parametrize("seq", ["AGTTNNN", "AGCTNNNNNNN", "AGCTNNNN", "AGCTNNN", "AGCTNNNN"])
    def test_gen_nearby_seqs_withn_none(self, maxdist, seq):
        """Test generation of nearby sequences."""

        # Setup
        barcode_sets = seq

        # Test
        nearby_seqs = list(gen_nearby_seqs(seq, barcode_sets[0], maxdist))

        # Assert
        assert len(nearby_seqs) == 0

    @pytest.mark.parametrize("maxdist", [1])
    @pytest.mark.parametrize(
        "seq, expected",
        [
            ("AGT", {"AAT", "ATT", "ACT", "TGT", "AGC", "AGA", "AGG", "CGT", "GGT"}),
            ("AGCT", {"TGCT", "GGCT", "CGCT", "ATCT", "AACT", "ACCT", "AGTT", "AGAT", "AGGT", "AGCC", "AGCA", "AGCG"}),
        ],
    )
    def test_gen_nearby_seqs_valid_hamm1(self, maxdist, seq, expected):
        """Test generation of nearby sequences."""

        # Setup
        barcode_sets = seq

        # Test
        nearby_seqs = set(gen_nearby_seqs(seq, barcode_sets[0], maxdist))

        # Assert
        assert_that(nearby_seqs).is_equal_to(expected)

    @pytest.mark.parametrize("maxdist", [2])
    @pytest.mark.parametrize(
        "seq, len_expected",
        [("AGT", 36), ("TNNCG", 16)],
    )
    def test_gen_nearby_seqs_valid_hamm2(self, maxdist, seq, len_expected):
        """Test generation of nearby sequences."""

        # Setup
        barcode_sets = seq

        # Test
        nearby_seqs = set(gen_nearby_seqs(seq, barcode_sets[0], maxdist))

        # Assert
        assert_that(len(nearby_seqs)).is_equal_to(len_expected)

    def test_generate_nearby_barcodes_by_length_none(self):
        assert_that(generate_nearby_barcodes_by_length).raises(TypeError).when_called_with(None, 1)

    def test_generate_nearby_barcodes_by_length_emptydict(self):
        assert_that(generate_nearby_barcodes_by_length({}, 1)).is_equal_to({})

    def test_generate_nearby_barcodes_by_length_invalid(self):
        # Call the function with a dictionary with invalid barcode values
        assert_that(generate_nearby_barcodes_by_length).raises(TypeError).when_called_with({"sample_1": 1234}, 1)

        # Call the function with an input that is not a dict
        assert_that(generate_nearby_barcodes_by_length).raises(TypeError).when_called_with("not_a_dict", 1)

    def test_generate_nearby_barcodes_by_length_valid(self):
        # Test with a valid barcode and distance
        grouped_barcode_dict = {
            4: {"sample_1": "ACGT"},
            8: {"sample_2": "ACG,AGG"},
        }
        distance = 1
        expected_result = {
            4: {"sample_1": {"ACCT", "ATGT", "ACGC", "ACGT", "ACAT", "CCGT", "GCGT", "ACGA", "ACGG", "AGGT", "TCGT", "AAGT", "ACTT"}},
            8: {
                "sample_2": {
                    "CCGAGG",
                    "ACGACG",
                    "ACGAGA",
                    "ACGAAG",
                    "TCGAGG",
                    "ACCAGG",
                    "ACGAGG",
                    "GCGAGG",
                    "ACAAGG",
                    "AGGAGG",
                    "ACGCGG",
                    "ACGAGT",
                    "ACGTGG",
                    "ACGATG",
                    "ACTAGG",
                    "ATGAGG",
                    "AAGAGG",
                    "ACGGGG",
                    "ACGAGC",
                }
            },
        }

        assert_that(generate_nearby_barcodes_by_length(grouped_barcode_dict, distance)).is_equal_to(expected_result)

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
        "fastq_file, barcode_sample_dict, max_hamming_distance, expected_samples, expected_file_content",
        [
            (
                "tests/data/seq/L002_R1.fastq",
                {
                    "sample_1": "ACTT,NTAT",
                    "sample_2": "ACGT,AGGT",
                },
                1,
                ["sample_1", "sample_2", "undetermined"],
                {
                    "sample_1": "GNAGGGGCGGCCCGGCCCCCACCCCCACGCCCGCCCGGGAGGCGGACGGGGGGAGAGGGAGAGCGCGGCGACGGGTATCTGGCTTCCTCGGCCCCGGGATTCGGCGAAAGCTGCGGCCGGAGGGCTGTAACACTCGGGGTGAGGTGGTAGA",
                    "sample_2": "",
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
            )
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
