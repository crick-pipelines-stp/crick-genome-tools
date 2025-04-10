"""
Tests for multiqc parser.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring

from assertpy import assert_that

from crick_genome_tools.reporting.multiqc.multiqc_parser import parse_samtools_stats


class TestMultiqc:
    def test_reporting_multiqc_parse_samtools_stats(self):
        # Setup and test
        df = parse_samtools_stats("tests/data/reporting/aav/samtools_host", [".host"])

        # Assert
        assert_that(df).is_not_none()
        assert_that(df.columns).contains("Sample")
