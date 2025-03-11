"""
Tests for toulligqc.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring


import unittest

from crick_genome_tools.reporting.multiqc.multiqc_parser import parse_samtools_stats


class TestMultiqc(unittest.TestCase):
    def test_multiqc_parse_samtools_stats(self):
        parse_samtools_stats("tests/data/reporting/aav/samtools_host", [".host"])
        raise NotImplementedError("Test not implemented")
