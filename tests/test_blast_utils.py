"""
Tests for fasta blast utils.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring

import os
from assertpy import assert_that
import pytest

from crick_genome_tools.blast.utils import build_fasta_from_top_hits

FASTA_SOURCE_FOLDER = "tests/data/io/fasta"
BLAST_FOLDER_VALID = "tests/data/blast/valid"
BLAST_FOLDER_INVALID = "tests/data/blast/invalid"

class TestBuildFastaFromTopHits:
    def test_blast_build_fasta_from_top_hits_empty_blast_file(self):
        # Test and Assert
        with pytest.raises(ValueError):
            assert_that(build_fasta_from_top_hits(os.path.join(FASTA_SOURCE_FOLDER, "extended_valid.fasta"), os.path.join(BLAST_FOLDER_INVALID, "empty.txt")))

    def test_blast_build_fasta_from_top_hits_blast_not_found(self):
        # Test and Assert
        with pytest.raises(FileNotFoundError):
            build_fasta_from_top_hits(os.path.join(FASTA_SOURCE_FOLDER, "extended_valid.fasta"), os.path.join(BLAST_FOLDER_INVALID, "not_found.txt"))

    def test_blast_build_fasta_from_top_hits_fasta_not_found(self):
        # Test and Assert
        with pytest.raises(FileNotFoundError):
            build_fasta_from_top_hits(os.path.join(FASTA_SOURCE_FOLDER, "not_found.fasta"), os.path.join(BLAST_FOLDER_VALID, "viral_top_hits.txt"))

    def test_blast_build_fasta_from_top_hits_single(self):
        # Test
        result = build_fasta_from_top_hits(
            os.path.join(FASTA_SOURCE_FOLDER, "extended_valid.fasta"), os.path.join(BLAST_FOLDER_VALID, "viral_top_hits_single.txt")
        )

        # Assert
        assert_that(result).is_not_empty()
        assert_that(result).contains("H1N1pdm_PA", "H1N1pdm_HA", "H1N1pdm_NA", "H1N1pdm_NS")
        assert_that(result).does_not_contain("fake_1", "fake_2")

    def test_blast_build_fasta_from_top_hits_multi(self):
        # Test
        result = build_fasta_from_top_hits(
            os.path.join(FASTA_SOURCE_FOLDER, "extended_valid.fasta"), os.path.join(BLAST_FOLDER_VALID, "viral_top_hits_multi.txt")
        )

        # Assert
        assert_that(result).is_not_empty()
        assert_that(result).contains("H1N1pdm_PA", "H1N1pdm_HA", "H1N1pdm_NA", "H1N1pdm_NS")

    def test_blast_build_fasta_from_top_hits_hit_not_found(self):
        # Test and Assert
        with pytest.raises(ValueError):
            build_fasta_from_top_hits(
                os.path.join(FASTA_SOURCE_FOLDER, "extended_valid.fasta"), os.path.join(BLAST_FOLDER_INVALID, "hit_not_found.txt")
            )
