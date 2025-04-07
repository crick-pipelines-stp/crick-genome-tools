"""
Tests covering the barcode_demux module
"""

# pylint: disable=missing-function-docstring,missing-class-docstring,no-member

from assertpy import assert_that

from crick_genome_tools.seq.barcode_demux import extract_index_from_header


class TestBarcodeDemux:
    def test_extract_index_from_header_none(self):
        assert_that(extract_index_from_header).raises(ValueError).when_called_with(None)

    def test_extract_index_from_header_empty(self):
        assert_that(extract_index_from_header).raises(ValueError).when_called_with("")

    def test_extract_index_from_header_invalid(self):
        header = "@LH00442:107:22YHM5LT3:2:1101:1092:1064 2:N:0:TCACCA:GGAC+NCCTTGT:CTC"
        assert_that(extract_index_from_header).raises(ValueError).when_called_with(header)

    def test_extract_index_from_header_isvalid(self):
        header = "@LH00442:107:22YHM5LT3:2:1101:1092:1064 2:N:0:TCACCAGGAC+NCCTTGTCTC"
        header_2 = "@LH00442:107 1:N:0:GCGCTTCTAC+NTCCTTGGCT"
        header_3 = "invalid_header"

        assert_that(extract_index_from_header(header)).is_equal_to("TCACCAGGAC+NCCTTGTCTC")
        assert_that(extract_index_from_header(header_2)).is_equal_to("GCGCTTCTAC+NTCCTTGGCT")
        assert_that(extract_index_from_header(header_3)).is_equal_to("")
