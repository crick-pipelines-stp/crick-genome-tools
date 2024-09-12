"""
Tests for iterative alignment
"""

# pylint: disable=missing-function-docstring,missing-class-docstring


import unittest

from crick_genome_tools.workflows.iterative_alignment import Aligner, IterationMode, IterativeAlignment
from tests.utils import with_temporary_folder


class TestLogSubprocess(unittest.TestCase):

    @with_temporary_folder
    def test_iterative_alignment(self, tmp_path):
        iter_align = IterativeAlignment(4, "/Users/cheshic/dev/test_data/iter_alignment", 1, 1, Aligner.BWA, IterationMode.COUNT, num_iterations=1)
        iter_align.run_sample(
            "sample_01",
            "tests/data/workflows/iter_align/fastq/random_10K/illumina_pe_nextera_test_10k_R1.fastq.gz",
            "tests/data/workflows/iter_align/fastq/random_10K/illumina_pe_nextera_test_10k_R2.fastq.gz",
            "tests/data/workflows/iter_align/refs/seg_ref_1.fasta",
        )
        raise NotImplementedError("Test not implemented")
