"""
Tests for toulligqc.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring


import unittest

from crick_genome_tools.reporting.tqc.configuration import ToulligqcConf
from crick_genome_tools.reporting.tqc.fastq_extractor import FastqExtractor


class TestToilligqc(unittest.TestCase):

    def test_toulligqc_extract_fastq_results_dict(self):
        # Setup
        config = ToulligqcConf()
        config["fastq"] = "tests/data/reporting/aav/toulligqc/P220.sub.1000.fastq.gz"
        config["images_directory"] = "data/reporting/aav/images"
        config["threshold"] = "10"
        config["batch_size"] = "1000"
        config["thread"] = "10"
        config["barcoding"] = "False"
        config["quiet"] = "False"
        config["barcode_selection"] = []

        # Test
        extractor = FastqExtractor(config)
        extractor.init()
        result_dict = {}
        extractor.extract(result_dict)

        # Assert
        self.assertIsInstance(result_dict, dict)
        self.assertGreater(len(result_dict), 0)
