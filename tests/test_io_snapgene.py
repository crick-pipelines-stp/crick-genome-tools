"""
Tests for using snapgene files.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring

from assertpy import assert_that

from crick_genome_tools.io.snapgene import convert_to_gff


class TestSnapGene:
    def test_io_snapgene_convert_to_gff_no_source(self, tmp_path):
        # Setup
        input_file = "tests/data/io/snapgene/P281_AAV.txt"
        output_file = str(tmp_path / "test.gff")
        contig_name = "contig_1"

        # Test
        convert_to_gff(input_file, output_file, contig_name)

        # Assert
        with open(output_file, "r", encoding="UTF-8") as f:
            lines = f.readlines()
            target_line = lines[0].split("\t")
            assert_that(target_line[0]).is_equal_to("contig_1")
            assert_that(target_line[1]).is_equal_to("SnapGene")
            assert_that(target_line[2]).is_equal_to("repeat_region")
            assert_that(target_line[3]).is_equal_to("1")
            assert_that(target_line[4]).is_equal_to("130")
            assert_that(target_line[5]).is_equal_to(".")
            assert_that(target_line[6]).is_equal_to("+")
            assert_that(target_line[7]).is_equal_to(".")
            # assert_that(target_line[8]).is_equal_to("Note=Functional equivalent of wild-type AAV2 ITR_AAV")

    def test_io_snapgene_convert_to_gff_with_source(self, tmp_path):
        # Setup
        input_file = "tests/data/io/snapgene/P281_AAV.txt"
        output_file = str(tmp_path / "test.gff")
        contig_name = "contig_1"
        source = "TestSource"

        # Test
        convert_to_gff(input_file, output_file, contig_name, source)

        # Assert
        with open(output_file, "r", encoding="UTF-8") as f:
            lines = f.readlines()
            target_line = lines[0].split("\t")
            assert_that(target_line[1]).is_equal_to("TestSource")
