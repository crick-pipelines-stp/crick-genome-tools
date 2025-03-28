"""
Tests reporting custom.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring

import os

from assertpy import assert_that

from crick_genome_tools.reporting.custom.samtools_parser import parse_samtools_flagstat, parse_samtools_idxstats
from crick_genome_tools.reporting.custom.mosdepth_parser import parse_mosdepth_per_base

class TestSamtoolsParser:
    def test_parse_samtools_flagstat(self, tmp_path):
        sample_clean = "sample_"


        # Create a temporary flagstat file
        flagstat_content = """123018 + 0 in total (QC-passed reads + QC-failed reads)
122215 + 0 primary
180 + 0 secondary
623 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
1482 + 0 mapped (1.20% : N/A)
679 + 0 primary mapped (0.56% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)"""
        flagstat_file = os.path.join(tmp_path, "sample_test.flagstat")
        with open(flagstat_file, "w", encoding="UTF-8") as f:
            f.write(flagstat_content)

        # Run the parser
        result = parse_samtools_flagstat(tmp_path, sample_clean)

        # Check the results
        assert_that(result).is_not_none()
        assert_that("test").is_in(result)
        assert_that(result["test"]["total_aligned"]).is_equal_to(123018)
        assert_that(result["test"]["primary"]).is_equal_to(122215)
        assert_that(result["test"]["secondary"]).is_equal_to(180)
        assert_that(result["test"]["supplementary"]).is_equal_to(623)
        assert_that(result["test"]["duplicates"]).is_equal_to(0)
        assert_that(result["test"]["primary_duplicates"]).is_equal_to(0)
        assert_that(result["test"]["mapped"]).is_equal_to(1482)
        assert_that(result["test"]["primary_mapped"]).is_equal_to(679)


    def test_parse_samtools_idxstats(self, tmp_path):
        sample_clean = "sample_"

        # Create a temporary flagstat file
        idxstat_content = """P281_AAV_5256bp	5256	98964	0
    pHelper	11635	18121	0
    aav2_1	7409	20969	0
    *	0	0	1257"""
        idxstat_file = os.path.join(tmp_path, "sample_test.idxstats")
        with open(idxstat_file, "w", encoding="UTF-8") as f:
            f.write(idxstat_content)

        # Run the parser
        result = parse_samtools_idxstats(tmp_path, sample_clean)

        # Check the results
        assert "test" in result
        assert result["test"]["P281_AAV_5256bp"]["length"] == 5256
        assert result["test"]["pHelper"]["reads"] == 18121

    def test_mosdepth_per_base(self):
        parse_mosdepth_per_base("tests/data/reporting/aav/coverage")
        raise NotImplementedError
