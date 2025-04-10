"""
This file contains unit tests for the Fasta class in the crick_genome_tools.io.fasta module.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring


import os
import pytest
from assertpy import assert_that

from crick_genome_tools.io.fasta import Fasta

TEST_FASTA_SOURCE_FOLDER = "tests/data/io/fasta"


class TestFasta:
    def test_fasta_read_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            Fasta.read_fasta_file(os.path.join(TEST_FASTA_SOURCE_FOLDER, "not_found.fasta"))

    def test_fasta_read_small_fasta(self):
        fasta_seqs = Fasta.read_fasta_file(os.path.join(TEST_FASTA_SOURCE_FOLDER, "H1N1pdm_HA.fasta"))
        assert_that(fasta_seqs).contains("H1N1pdm_HA")
        assert_that(fasta_seqs["H1N1pdm_HA"]).is_length(1701)

    def test_fasta_read_fasta_is_uppercase(self):
        fasta_seqs = Fasta.read_fasta_file(os.path.join(TEST_FASTA_SOURCE_FOLDER, "H1N1pdm_HA.fasta"))
        assert_that(fasta_seqs["H1N1pdm_HA"]).is_equal_to(fasta_seqs["H1N1pdm_HA"].upper())

    def test_fasta_read_fasta_with_missing_tag_raises(self):
        with pytest.raises(ValueError):
            Fasta.read_fasta_file(os.path.join(TEST_FASTA_SOURCE_FOLDER, "invalid", "no_tag.fasta"))

    def test_fasta_read_directory_does_not_exist(self):
        with pytest.raises(FileNotFoundError):
            Fasta.read_fasta_directory(os.path.join(TEST_FASTA_SOURCE_FOLDER, "test"))

    def test_fasta_read_directory_success(self):
        fasta_seqs = Fasta.read_fasta_directory(os.path.join(TEST_FASTA_SOURCE_FOLDER, "valid_dir"))
        assert_that(fasta_seqs).is_length(4)

    def test_fasta_read_multi_tag(self):
        fasta_seqs = Fasta.read_fasta_file(os.path.join(TEST_FASTA_SOURCE_FOLDER, "multi_tag.fasta"))
        assert_that(fasta_seqs).is_length(2)

    def test_fasta_write_fasta_to_disk(self, tmp_path):
        fasta_seqs = Fasta.read_fasta_file(os.path.join(TEST_FASTA_SOURCE_FOLDER, "H1N1pdm_HA.fasta"))
        output_file = tmp_path / "H1N1pdm_HA.fasta"
        Fasta.write_fasta_file(fasta_seqs, str(output_file))

        assert_that(output_file.exists()).is_true()
        assert_that(output_file.read_text("utf-8")).contains("H1N1pdm_HA")

    @pytest.mark.parametrize(
        "ref_file, key, seq",
        [
            ("short_lower", "lower", "ACCT"),
            ("short_upper", "upper", "ACCT"),
        ]
    )
    def test_fasta_read_parametrized_sequences(self, ref_file, key, seq):
        fasta_seqs = Fasta.read_fasta_file(os.path.join(TEST_FASTA_SOURCE_FOLDER, f"{ref_file}.fasta"))
        assert_that(fasta_seqs).contains(key)
        assert_that(fasta_seqs[key]).is_equal_to(seq)
