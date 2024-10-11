"""
Tests for iterative alignment
"""

# pylint: disable=missing-function-docstring,missing-class-docstring


import os
import unittest

import pytest

from crick_genome_tools.io.fasta import Fasta
from crick_genome_tools.workflows.iterative_alignment import Aligner, IterationMode, IterativeAlignment


class TestIterativeAlignment(unittest.TestCase):

    @pytest.mark.container(image="thecrick/pipetech_iterative_alignment:test", test_dir="/Users/cheshic/dev/test_data/iter_alignment")
    def test_iterative_alignment(self):
        if not os.getenv("IN_CONTAINER"):
            return

        # Setup
        # output_path = os.environ.get("TMPDIR")
        output_path = "/Users/cheshic/dev/test_data/iter_alignment"
        iter_align = IterativeAlignment(
            num_cores=8,
            output_path=output_path,
            min_iterations=1,
            max_iterations=2,
            aligner=Aligner.BWA,
            iteration_mode=IterationMode.COUNT,
            var_thresh=0.001,
            min_coverage=100,
            bwa_args=["-T10"],
        )
        sample_id = "sample_01"
        ref_file = "tests/data/workflows/iter_align/refs/seg_ref_1.fasta"
        r1_file = "tests/data/workflows/iter_align/fastq/random_10K/illumina_pe_nextera_test_10k_R1.fastq.gz"
        r2_file = "tests/data/workflows/iter_align/fastq/random_10K/illumina_pe_nextera_test_10k_R2.fastq.gz"

        # Test
        iter_align.run_sample(
            sample_id,
            r1_file,
            r2_file,
            ref_file,
        )
        # raise NotImplementedError("MANUAL END")

        # Verify that the final consensus file was created
        final_consensus_path = os.path.join(output_path, sample_id, f"{sample_id}_final_consensus.fasta")
        assert os.path.exists(final_consensus_path), "Final consensus FASTA not found."

        # Read and verify that the contents of the final consensus file is the same as the reference
        # final_consensus = Fasta.read_fasta_file(final_consensus_path)
        # ref_seq = Fasta.read_fasta_file(ref_file)
        # assert final_consensus == ref_seq, "Final consensus does not match reference."

        # Verify the BAM and flagstat files
        final_bam_path = os.path.join(output_path, sample_id, f"{sample_id}_final.bam")
        assert os.path.exists(final_bam_path), "Final BAM file not found."

        final_flagstat_path = os.path.join(output_path, sample_id, f"{sample_id}_final.flagstat")
        assert os.path.exists(final_flagstat_path), "Final flagstat file not found."

        # Verify outputs of each iteration
        for i in range(1, 3):
            iteration_dir = os.path.join(output_path, sample_id, f"iteration_{i}")
            assert os.path.exists(iteration_dir), f"Iteration directory {iteration_dir} not found."

            # Check for consensus FASTA in each iteration
            consensus_fasta_path = os.path.join(iteration_dir, f"{sample_id}_iter_{i}.consensus.fasta")
            assert os.path.exists(consensus_fasta_path), f"Consensus FASTA for iteration {i} not found."

            # Read and verify the consensus FASTA
            iter_consensus_fasta = Fasta.read_fasta_file(consensus_fasta_path)
            assert len(iter_consensus_fasta) > 0, f"Consensus FASTA for iteration {i} is empty."

            # Check Hamming distance log
            log_file = os.path.join(iteration_dir, "logs", f"{sample_id}_iter_{i}.consensus.log")
            assert os.path.exists(log_file), f"Log file for iteration {i} not found."
