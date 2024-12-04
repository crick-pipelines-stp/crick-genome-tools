"""
Tests for samtools utils.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring

import os
import unittest
from tests.utils import with_temporary_folder

from crick_genome_tools.samtools.utils import count_table_from_pileup

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
