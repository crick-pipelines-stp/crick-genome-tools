"""
Tests reporting custom.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring

from crick_genome_tools.reporting.custom.samtools_parser import parse_samtools_flagstat
import os

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
        assert "test" in result
        assert result["test"]["total_aligned"] == 123018
        assert result["test"]["primary"] == 122215
        assert result["test"]["secondary"] == 180
        assert result["test"]["supplementary"] == 623
        assert result["test"]["duplicates"] == 0
        assert result["test"]["primary_duplicates"] == 0
        assert result["test"]["mapped"] == 1482
        assert result["test"]["primary_mapped"] == 679
