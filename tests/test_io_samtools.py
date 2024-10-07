"""
Tests for samtools.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring


from io import StringIO

import pytest

from crick_genome_tools.io.samtools import parse_flagstat


class TestSamtools:

    # Mock function to simulate reading from a file
    def mock_open_flagstat(self, mock_content):
        return StringIO(mock_content)

    def test_parse_flagstat_normal_case(self, monkeypatch):
        # Simulate a typical flagstat file
        mock_flagstat_content = """100000 + 0 in total (QC-passed reads + QC-failed reads)
    50000 + 0 secondary
    0 + 0 supplementary
    0 + 0 duplicates
    80000 + 0 mapped (80.00%)
    """
        monkeypatch.setattr("builtins.open", lambda *args, **kwargs: self.mock_open_flagstat(mock_flagstat_content))

        total_reads, mapped_reads, alignment_rate = parse_flagstat("fake_file.txt")
        assert total_reads == 100000
        assert mapped_reads == 80000
        assert alignment_rate == 80.00

    def test_parse_flagstat_zero_reads_case(self, monkeypatch):
        # Simulate a flagstat file with zero total reads
        mock_flagstat_content = """0 + 0 in total (QC-passed reads + QC-failed reads)
    0 + 0 secondary
    0 + 0 supplementary
    0 + 0 duplicates
    0 + 0 mapped (0.00%)
    """
        monkeypatch.setattr("builtins.open", lambda *args, **kwargs: self.mock_open_flagstat(mock_flagstat_content))

        total_reads, mapped_reads, alignment_rate = parse_flagstat("fake_file.txt")
        assert total_reads == 0
        assert mapped_reads == 0
        assert alignment_rate == 0.00

    def test_parse_flagstat_incorrect_format_case(self, monkeypatch):
        # Simulate a flagstat file with an unexpected format
        mock_flagstat_content = """Some unexpected content
    Another random line
    0 + 0 supplementary
    Not the right line
    0 + 0 mapped (100.00%)
    """
        monkeypatch.setattr("builtins.open", lambda *args, **kwargs: self.mock_open_flagstat(mock_flagstat_content))

        with pytest.raises(ValueError):
            parse_flagstat("fake_file.txt")
