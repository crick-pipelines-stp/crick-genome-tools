"""
Tests for report data parser.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring

import os
from assertpy import assert_that
from crick_genome_tools.reporting.report_data_parser import ReportDataParser

class TestReportDataParser:

    def test_rdp_init(self):
        parser = ReportDataParser("tests/data/reporting/aav")
        assert_that(parser.data_folder).is_equal_to("tests/data/reporting/aav")
        assert_that(parser.result_dict).is_instance_of(dict)
        assert_that(parser.dataframe_dict).is_instance_of(dict)
        assert_that(parser.folder_names).is_equal_to(os.listdir("tests/data/reporting/aav"))

    def test_rdp_get_data_aav(self):
        # Setup
        parser = ReportDataParser("tests/data/reporting/aav")

        # Test
        parser.get_data()

        # Assert
        assert_that(len(parser.result_dict)).is_equal_to(3)
        assert_that(len(parser.dataframe_dict)).is_equal_to(3)

    def test_rdp_get_samtools_host_data(self):
        # Setup
        parser = ReportDataParser("tests/data/reporting/aav")

        # Test
        parser.get_samtools_host_data("tests/data/reporting/aav/samtools_host", ".host", "host")

        # Assert
        assert_that(parser.merged_dataframe_dict.keys()).contains("samtools_host")

    def test_rdp_get_samtools_contam_data(self):
        # Setup
        parser = ReportDataParser("tests/data/reporting/aav")
        parser.get_data()

        # Test
        # parser.get_samtools_contam_data("tests/data/reporting/aav/samtools_contaminent", ".contam", "contam")

        # Assert
        # assert_that(parser.merged_dataframe_dict.keys()).contains("samtools_host")
        raise NotImplementedError