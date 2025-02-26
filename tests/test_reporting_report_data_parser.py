"""
Tests for report data parser.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring

import unittest
import os
from crick_genome_tools.reporting.report_data_parser import ReportDataParser

class TestReportDataParser(unittest.TestCase):

    def test_rdp_init(self):
        parser = ReportDataParser("tests/data/reporting/aav")
        self.assertEqual(parser.data_folder, "tests/data/reporting/aav")
        self.assertIsInstance(parser.result_dict, dict)
        self.assertIsInstance(parser.dataframe_dict, dict)
        self.assertEqual(parser.folder_names, os.listdir("tests/data/reporting/aav"))

    def test_rdp_get_data_aav(self):
        parser = ReportDataParser("tests/data/reporting/aav")
        parser.get_data()

        self.assertEqual(len(parser.result_dict), 3)
        self.assertEqual(len(parser.dataframe_dict), 3)
