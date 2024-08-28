"""
Tests for fasta file reading
"""

# pylint: disable=missing-function-docstring,missing-class-docstring

import os
import unittest

import pytest

from crick_genome_tools.io.fasta import Fasta
from tests.utils import with_temporary_folder


TEST_FASTA_SOURCE_FOLDER = "tests/data/io/fasta"


class TestIoFasta(unittest.TestCase):
    def test_fasta_read_file_not_found(self):
        # Test and Assert
        with self.assertRaises(FileNotFoundError):
            Fasta.read_fasta_file(os.path.join(TEST_FASTA_SOURCE_FOLDER, "not_found.fasta"))

    def test_fasta_read_small_fasta(self):
        # Test
        fasta_seqs = Fasta.read_fasta_file(os.path.join(TEST_FASTA_SOURCE_FOLDER, "H1N1pdm_HA.fasta"))

        # Assert
        self.assertTrue("H1N1pdm_HA" in fasta_seqs)
        self.assertTrue(len(fasta_seqs["H1N1pdm_HA"]) == 1701)

    def test_fasta_read_fasta_upper(self):
        # Test
        fasta_seqs = Fasta.read_fasta_file(os.path.join(TEST_FASTA_SOURCE_FOLDER, "H1N1pdm_HA.fasta"))

        # Assert
        self.assertEqual(fasta_seqs["H1N1pdm_HA"], fasta_seqs["H1N1pdm_HA"].upper())

    def test_fasta_read_notag_error(self):
        # Test and Assert
        with self.assertRaises(ValueError):
            Fasta.read_fasta_file(os.path.join(TEST_FASTA_SOURCE_FOLDER, "invalid", "no_tag.fasta"))

    def test_fasta_read_dir_doesnotexist(self):
        # Test and Assert
        with self.assertRaises(FileNotFoundError):
            Fasta.read_fasta_directory(os.path.join(TEST_FASTA_SOURCE_FOLDER, "test"))

    def test_fasta_read_dir(self):
        # Test
        fasta_seqs = Fasta.read_fasta_directory(os.path.join(TEST_FASTA_SOURCE_FOLDER, "valid_dir"))

        # Assert
        self.assertTrue(len(fasta_seqs) == 4)

    def test_fasta_multi_tag(self):
        # Test
        fasta_seqs = Fasta.read_fasta_file(os.path.join(TEST_FASTA_SOURCE_FOLDER, "multi_tag.fasta"))

        # Assert
        self.assertTrue(len(fasta_seqs) == 2)

    @with_temporary_folder
    def test_fasta_write(self, tmpdirname):
        # Test
        fasta_seqs = Fasta.read_fasta_file(os.path.join(TEST_FASTA_SOURCE_FOLDER, "H1N1pdm_HA.fasta"))
        Fasta.write_fasta_file(fasta_seqs, os.path.join(tmpdirname, "H1N1pdm_HA.fasta"))

        # Assert
        self.assertTrue(os.path.exists(os.path.join(tmpdirname, "H1N1pdm_HA.fasta")))
        with open(os.path.join(tmpdirname, "H1N1pdm_HA.fasta"), "r", encoding="UTF-8") as file:
            self.assertEqual(file.read(), ">H1N1pdm_HA\n" + fasta_seqs["H1N1pdm_HA"] + "\n")


class TestFastaFixture:
    @pytest.mark.parametrize("ref_file, key, seq", [("short_lower", "lower", "ACCT"), ("short_upper", "upper", "ACCT")])
    def test_fasta_read_with_expected(self, ref_file, key, seq):  # pylint: disable=W0613
        # Test
        fasta_seqs = Fasta.read_fasta_file(os.path.join(TEST_FASTA_SOURCE_FOLDER, ref_file + ".fasta"))

        # Assert
        assert key in fasta_seqs
        assert fasta_seqs[key] == seq
