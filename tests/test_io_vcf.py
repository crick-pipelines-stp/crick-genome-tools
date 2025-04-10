"""
Tests vcf utils
"""

# pylint: disable=missing-function-docstring,missing-class-docstring


import pytest
from assertpy import assert_that

from crick_genome_tools.io.vcf import determine_variant_type, generate_merged_vcf_report


class TestVcf:
    def test_io_vcf_determine_variant_type_snp(self):
        assert_that(determine_variant_type("A", "T")).is_equal_to("SNV")

    def test_io_vcf_determine_variant_type_ins(self):
        assert_that(determine_variant_type("A", "AT")).is_equal_to("INS")

    def test_io_vcf_determine_variant_type_del(self):
        assert_that(determine_variant_type("AT", "A")).is_equal_to("DEL")

    def test_io_vcf_determine_variant_type_indel(self):
        assert_that(determine_variant_type("AT", "GC")).is_equal_to("INDEL")

    def test_generate_merged_vcf_report_bad_vcf(self):
        vcf_files = [
            "tests/data/io/vcf/error.vcf",
            "tests/data/io/vcf/FAY66992_BC15.pass.vcf",
        ]
        with pytest.raises(ValueError):
            generate_merged_vcf_report(vcf_files, ["test"])

    def test_io_vcf_generate_merged_vcf_report_insufficient_files(self):
        vcf_files = ["tests/data/io/vcf/FAY66992_BC15.pass.vcf"]
        with pytest.raises(ValueError):
            generate_merged_vcf_report(vcf_files, ["medaka"])

    def test_io_vcf_generate_merged_vcf_report_tool_mismatch(self):
        vcf_files = [
            "tests/data/io/vcf/FAY66992_BC15.pass.vcf",
            "tests/data/io/vcf/FAY66992_BC15.clair3.merge_output.vcf",
        ]
        tool_names = ["medaka"]
        with pytest.raises(ValueError):
            generate_merged_vcf_report(vcf_files, tool_names)

    def test_io_vcf_generate_merged_vcf_report_nanopore(self):
        vcf_files = [
            "tests/data/io/vcf/FAY66992_BC15.clair3.merge_output.vcf",
            "tests/data/io/vcf/FAY66992_BC15.pass.vcf",
            "tests/data/io/vcf/FAY66992_BC15.lofreq.vcf",
            "tests/data/io/vcf/FAY66992_BC15.sniffles.vcf",
            "tests/data/io/vcf/FAY66992_BC15.snpeff.vcf",
        ]
        variants, _, _ = generate_merged_vcf_report(vcf_files, ["clair3", "medaka", "lofreq", "sniffles", "snpeff"])
        assert_that(variants).is_length(148)

    def test_io_vcf_generate_merged_vcf_report_illumina(self):
        vcf_files = [
            "tests/data/io/vcf/20-A_Tajikistan_12-928_2023.lofreq.vcf",
            "tests/data/io/vcf/20-A_Tajikistan_12-928_2023.freebayes.vcf",
            "tests/data/io/vcf/20-A_Tajikistan_12-928_2023.snpeff.vcf",
        ]
        variants, _, _ = generate_merged_vcf_report(vcf_files, ["lofreq", "freebayes", "snpeff"])
        assert_that(variants).is_length(260)

    # def test_generate_merged_vcf_report_dev(self):
    #     # Setup
    #     vcf_files = [
    #                  "tests/data/io/vcf/FAY66992_BC15.clair3.merge_output.vcf",
    #                  "tests/data/io/vcf/FAY66992_BC15.pass.vcf",
    #                  "tests/data/io/vcf/FAY66992_BC15.lofreq.vcf",
    #                  "tests/data/io/vcf/FAY66992_BC15.sniffles.vcf",
    #                  "tests/data/io/vcf/FAY66992_BC15.snpeff.vcf"
    #                 ]

    #     # Test
    #     variants, header, report = generate_merged_vcf_report(vcf_files, ["clair3", "medaka", "lofreq", "sniffles", "snpeff"], "output.tsv")

    #     raise NotImplementedError("Test not implemented")
