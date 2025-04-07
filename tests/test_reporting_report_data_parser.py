"""
Tests for report data parser.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring

import os

from assertpy import assert_that

from crick_genome_tools.reporting.report_data_parser import ReportDataParser


class TestReportDataParser:

    def test_rdp_init(self):
        # Init
        data_path = "tests/data/reporting/aav"
        parser = ReportDataParser(data_path)

        # Assert
        assert_that(parser.data_folder).is_equal_to(data_path)
        assert_that(parser.result_dict).is_instance_of(dict)
        assert_that(parser.dataframe_dict).is_instance_of(dict)
        assert_that(parser.folder_names).is_equal_to(os.listdir(data_path))

    def test_rdp_save_data(self, tmp_path):
        # Setup
        parser = ReportDataParser("tests/data/reporting/aav")
        output_file = os.path.join(tmp_path, "test.pkl")

        # Test
        parser.save_data(output_file)

        # Assert
        assert_that(os.path.exists(output_file)).is_true()
        # parser.save_data("tests/data/reporting/aav.pkl")

    def test_rdp_get_data_aav(self):
        # Setup
        parser = ReportDataParser("tests/data/reporting/aav")

        # Test
        parser.get_data()

        # Assert
        assert_that(len(parser.result_dict)).is_equal_to(3)
        assert_that(len(parser.dataframe_dict)).is_equal_to(3)

    def test_rdp_get_samtools_host_data(self):
        # Setup
        parser = ReportDataParser("tests/data/reporting/aav")

        # Test
        parser.get_samtools_host_data("tests/data/reporting/aav/samtools_host", ".host", "host")

        # Assert
        assert_that(parser.merged_dataframe_dict.keys()).contains("samtools_host")

    def test_rdp_get_samtools_contam_data(self):
        # Setup
        parser = ReportDataParser("tests/data/reporting/aav")

        # Test
        parser.get_samtools_contam_data("tests/data/reporting/aav/samtools_contaminent", ".contam", "contam")

        # Assert
        assert_that(parser.merged_dataframe_dict.keys()).contains("samtools_contam")
