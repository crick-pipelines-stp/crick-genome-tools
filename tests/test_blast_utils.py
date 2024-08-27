"""
Tests for fasta blast utils.
"""
# pylint: disable=missing-function-docstring,missing-class-docstring

import os
import unittest

from crick_genome_tools.blast.utils import build_fasta_from_top_hits

FASTA_SOURCE_FOLDER = "tests/data/io/fasta"
BLAST_FOLDER_VALID = "tests/data/blast/valid"
BLAST_FOLDER_INVALID = "tests/data/blast/invalid"

class TestBuildFastaFromTopHits(unittest.TestCase):
    def test_build_fasta_from_top_hits_empty_blast_file(self):
        # Test and Assert
        with self.assertRaises(ValueError):
            build_fasta_from_top_hits(os.path.join(FASTA_SOURCE_FOLDER, "extended_valid.fasta"),
                                      os.path.join(BLAST_FOLDER_INVALID, "empty.txt"))

    def test_build_fasta_from_top_hits_blast_not_found(self):
        # Test and Assert
        with self.assertRaises(FileNotFoundError):
            build_fasta_from_top_hits(os.path.join(FASTA_SOURCE_FOLDER, "extended_valid.fasta"),
                                      os.path.join(BLAST_FOLDER_INVALID, "not_found.txt"))

    def test_build_fasta_from_top_hits_fasta_not_found(self):
        # Test and Assert
        with self.assertRaises(FileNotFoundError):
            build_fasta_from_top_hits(os.path.join(FASTA_SOURCE_FOLDER, "not_found.fasta"),
                                      os.path.join(BLAST_FOLDER_VALID, "viral_top_hits.txt"))

    def test_build_fasta_from_top_hits_success(self):
        # Test
        result = build_fasta_from_top_hits(os.path.join(FASTA_SOURCE_FOLDER, "extended_valid.fasta"),
                                           os.path.join(BLAST_FOLDER_VALID, "viral_top_hits.txt"))

        # Assert
        self.assertTrue("H1N1pdm_PA" in result)
        self.assertTrue("H1N1pdm_HA" in result)
        self.assertTrue("H1N1pdm_NA" in result)
        self.assertTrue("H1N1pdm_NS" in result)
        self.assertTrue("fake_1" not in result)
        self.assertTrue("fake_2" not in result)

    def test_build_fasta_from_top_hits_duplicate_hits(self):
        # Test and Assert
        with self.assertRaises(ValueError):
            build_fasta_from_top_hits(os.path.join(FASTA_SOURCE_FOLDER, "extended_valid.fasta"),
                                      os.path.join(BLAST_FOLDER_INVALID, "duplicate_hits.txt"))

    def test_build_fasta_from_top_hits_hit_not_found(self):
        # Test and Assert
        with self.assertRaises(ValueError):
            build_fasta_from_top_hits(os.path.join(FASTA_SOURCE_FOLDER, "extended_valid.fasta"),
                                      os.path.join(BLAST_FOLDER_INVALID, "hit_not_found.txt"))
