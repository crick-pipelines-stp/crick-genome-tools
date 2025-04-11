"""
Tests covering the barcode_demux module
"""

# pylint: disable=missing-function-docstring,missing-class-docstring,no-member

from assertpy import assert_that

from crick_genome_tools.seq.barcode_demux import extract_index_from_header_illumina, generate_nearby_barcodes_by_length, group_samples_by_index_length


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

    def test_generate_nearby_barcodes_by_length_none(self):
        assert_that(generate_nearby_barcodes_by_length).raises(TypeError).when_called_with(None, 1)

    def test_generate_nearby_barcodes_by_length_emptydict(self):
        assert_that(generate_nearby_barcodes_by_length({}, 1)).is_equal_to([])

    def test_generate_nearby_barcodes_by_length_invalid(self):
        # Call the function with a dictionary with invalid barcode values
        assert_that(generate_nearby_barcodes_by_length).raises(TypeError).when_called_with({"sample_1": 1234}, 1)

        # Call the function with an input that is not a dict
        assert_that(generate_nearby_barcodes_by_length).raises(TypeError).when_called_with("not_a_dict", 1)

    def test_generate_nearby_barcodes_by_length(self):
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
                    "ACG,AGG",
                    "ACG,AGT",
                    "ACG,AGC",
                    "ACG,ACC",
                    "ACG,ATA",
                    "ACG,ATG",
                    "ACG,ATC",
                    "ACG,ATT",
                    "ACG,ACA",
                    "ACG,AAG",
                    "ACG,AAT",
                    "ACG,AAC",
                }
            },
        }

        assert_that(generate_nearby_barcodes_by_length(grouped_barcode_dict, distance)).is_equal_to(expected_result)
