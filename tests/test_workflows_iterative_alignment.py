"""
Tests for iterative alignment
"""

# pylint: disable=missing-function-docstring,missing-class-docstring


import os
import unittest

import pytest

from crick_genome_tools.workflows.iterative_alignment import Aligner, IterationMode, IterativeAlignment


class TestIterativeAlignment(unittest.TestCase):

    @pytest.mark.container(image="thecrick/pipetech_iterative_alignment:test", test_dir="/Users/cheshic/dev/test_data/iter_alignment")
    def test_iterative_alignment(self):
        iter_align = IterativeAlignment(
            num_cores=8,
            # output_path=os.environ.get('TMPDIR'),
            output_path="/Users/cheshic/dev/test_data/iter_alignment",
            min_iterations=1,
            max_iterations=2,
            aligner=Aligner.BWA,
            iteration_mode=IterationMode.COUNT,
            bwa_args=["-T10"],
        )

        iter_align.run_sample(
            "sample_01",
            "tests/data/workflows/iter_align/fastq/random_10K/illumina_pe_nextera_test_10k_R1.fastq.gz",
            "tests/data/workflows/iter_align/fastq/random_10K/illumina_pe_nextera_test_10k_R2.fastq.gz",
            "tests/data/workflows/iter_align/refs/seg_ref_1.fasta",
        )
        raise NotImplementedError("MANUAL END")
