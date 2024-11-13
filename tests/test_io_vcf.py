"""
Tests vcf utils
"""

# pylint: disable=missing-function-docstring,missing-class-docstring


import unittest

from crick_genome_tools.io.vcf import generate_merged_vcf_report
from crick_genome_tools.io.vcf import determine_variant_type


class TestVcf(unittest.TestCase):
    def test_determine_variant_type_snp(self):
        self.assertEqual(determine_variant_type("A", "T"), "SNV")


    def test_determine_variant_type_ins(self):
        self.assertEqual(determine_variant_type("A", "AT"), "INS")


    def test_determine_variant_type_del(self):
        self.assertEqual(determine_variant_type("AT", "A"), "DEL")
 

    def test_determine_variant_type_indel(self):
        self.assertEqual(determine_variant_type("AT", "GC"), "INDEL")


    def test_generate_merged_vcf_report_bad_vcf(self):
        # Setup
        vcf_files = ["tests/data/io/vcf/error.vcf", "tests/data/io/vcf/FAY66992_BC15.pass.vcf"]

        # Test and Assert
        with self.assertRaises(ValueError):
            generate_merged_vcf_report(vcf_files, ["test"])

    def test_generate_merged_vcf_report_insufficient_vcf_files(self):
        # Setup
        vcf_files = ["tests/data/io/vcf/FAY66992_BC15.pass.vcf"]

        # Test and Assert
        with self.assertRaises(ValueError):
            generate_merged_vcf_report(vcf_files, ["medaka"])


    def test_generate_merged_vcf_report_tool_number_mismatch(self):
        # Setup
        vcf_files = [
            "tests/data/io/vcf/FAY66992_BC15.pass.vcf",
            "tests/data/io/vcf/FAY66992_BC15.clair3.merge_output.vcf"
        ]
        tool_names = ["medaka"]

        # Test and Assert
        with self.assertRaises(ValueError):
            generate_merged_vcf_report(vcf_files, tool_names)


    def test_generate_merged_vcf_report(self):
        # Setup
        vcf_files = [
                     "tests/data/io/vcf/FAY66992_BC15.clair3.merge_output.vcf",
                     "tests/data/io/vcf/FAY66992_BC15.pass.vcf",
                     "tests/data/io/vcf/FAY66992_BC15.lofreq.vcf",
                     "tests/data/io/vcf/FAY66992_BC15.sniffles.vcf",
                     "tests/data/io/vcf/FAY66992_BC15.snpeff.vcf"
                    ]

        # Test
        variants, header, report = generate_merged_vcf_report(vcf_files, ["clair3", "medaka", "lofreq", "sniffles", "snpeff"], 20)

        # Assert variant count
        self.assertEqual(len(variants), 148)


    # def test_generate_merged_vcf_report_dev(self):
    #     # Setup
    #     vcf_files = [
    #                  "tests/data/io/vcf/FAY66992_BC15.pass.vcf",
    #                  "tests/data/io/vcf/FAY66992_BC15.clair3.merge_output.vcf",
    #                  "tests/data/io/vcf/FAY66992_BC15.lofreq.vcf",
    #                  "tests/data/io/vcf/FAY66992_BC15.snpeff.vcf",
    #                 ]

    #     # Test
    #     variants, header, report = generate_merged_vcf_report(vcf_files, ["medaka", "clair3", "lofreq", "snpeff"], 20, "output.tsv")

    #     raise NotImplementedError("Test not implemented")
