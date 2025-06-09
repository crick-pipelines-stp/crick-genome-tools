"""
Tests for samtools utils.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring

from assertpy import assert_that

from crick_genome_tools.samtools.utils import count_table_from_pileup


class TestSamtoolsUtils:
    def test_samtools_utils_count_table_from_pileup_empty_file(self, tmp_path):
        output_path = tmp_path / "count_table.txt"
        empty_input = tmp_path / "empty.txt"
        empty_input.write_text("", encoding="utf-8")

        count_table_from_pileup(str(empty_input), str(output_path))

        assert_that(output_path.exists()).is_true()
        lines = output_path.read_text(encoding="utf-8").splitlines()
        assert_that(lines).is_length(1)

    def test_samtools_utils_count_table_from_pileup_divide_by_zero(self, tmp_path):
        output_path = tmp_path / "count_table.txt"
        input_path = tmp_path / "divide_by_zero.txt"
        input_path.write_text("chr1\t1\tN\t0\t*\t*\n", encoding="utf-8")

        count_table_from_pileup(str(input_path), str(output_path))

        assert_that(output_path.exists()).is_true()
        lines = output_path.read_text(encoding="utf-8").splitlines()
        assert_that(lines).is_length(2)

    def test_samtools_utils_count_table_from_pileup_valid_single_contig(self, tmp_path):
        output_path = tmp_path / "count_table.txt"

        count_table_from_pileup("tests/data/mpileup/h1n1.txt", str(output_path))

        assert_that(output_path.exists()).is_true()

        lines = output_path.read_text(encoding="utf-8").splitlines()
        assert_that(lines).is_length(401)

    def test_samtools_utils_count_table_from_pileup_valid_multi_contig(self, tmp_path):
        output_path = tmp_path / "count_table.tsv"

        count_table_from_pileup("tests/data/mpileup/SRR33681751_H5N1_cattle.mpileup", str(output_path))

        assert_that(output_path.exists()).is_true()

        lines = output_path.read_text(encoding="utf-8").splitlines()
        assert_that(lines).is_length(1501)
