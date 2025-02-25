"""
Tests for samtools utils.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring

import os
import unittest

from crick_genome_tools.samtools.utils import count_table_from_pileup
from tests.utils import with_temporary_folder


class TestSamtoolsUtils(unittest.TestCase):

    @with_temporary_folder
    def test_count_table_from_pileup(self, tmp_dir):
        # Setup
        output_path = os.path.join(tmp_dir, "count_table.txt")

        # Test
        count_table_from_pileup("tests/data/mpileup/h1n1.txt", output_path)

        # Assert
        self.assertTrue(os.path.exists(output_path))
        with open(output_path, "r", encoding="UTF-8") as output_file:
            output_lines = output_file.readlines()
            self.assertEqual(len(output_lines), 401)

    @with_temporary_folder
    def test_count_table_from_pileup_empty_file(self, tmp_dir):
        # Setup
        output_path = os.path.join(tmp_dir, "count_table.txt")
        empty_input_path = os.path.join(tmp_dir, "empty.txt")
        with open(empty_input_path, "w", encoding="UTF-8") as empty_file:
            empty_file.write("")

        # Test
        count_table_from_pileup(empty_input_path, output_path)

        # Assert
        self.assertTrue(os.path.exists(output_path))
        with open(output_path, "r", encoding="UTF-8") as output_file:
            output_lines = output_file.readlines()
            self.assertEqual(len(output_lines), 1)

    @with_temporary_folder
    def test_count_table_from_pileup_divide_by_zero(self, tmp_dir):
        # Setup
        output_path = os.path.join(tmp_dir, "count_table.txt")
        divide_by_zero_input_path = os.path.join(tmp_dir, "divide_by_zero.txt")
        with open(divide_by_zero_input_path, "w", encoding="UTF-8") as input_file:
            input_file.write("chr1\t1\tN\t0\t*\t*\n")

        # Test
        count_table_from_pileup(divide_by_zero_input_path, output_path)

        # Assert
        self.assertTrue(os.path.exists(output_path))
        with open(output_path, "r", encoding="UTF-8") as output_file:
            output_lines = output_file.readlines()
            self.assertEqual(len(output_lines), 2)
