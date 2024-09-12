"""
Performs iterative alignment of sequences to a reference genome.
"""

import logging
import os
import shutil
from enum import Enum
import subprocess

from crick_genome_tools.io.log_subprocess import LogSubprocess

log = logging.getLogger(__name__)

class IterationMode(Enum):
    """
    Enum for aligner to use
    """
    COUNT = "count"
    HAMMING = "hamming"


class Aligner(Enum):
    """
    Enum for aligner to use
    """
    BWA = "bwa"
    # BOWTIE2 = "bowtie2"
    # HISAT2 = "hisat2"


class IterativeAlignment:
    """
    Perform iterative alignment of sequences to a reference genome.
    """

    def __init__(self, num_cores: int, output_path: str, min_iterations: int, max_iterations: int, aligner: Aligner, iteration_mode: IterationMode, **kwargs):
        """
        Initialise the IterativeAlignment object
        """
        self.num_cores = num_cores
        self.output_path = output_path
        self.min_iterations = min_iterations
        self.max_iterations = max_iterations
        self.aligner = aligner
        self.iteration_mode = iteration_mode

        # Mode-specific configurations
        if iteration_mode == IterationMode.COUNT:
            self.num_iterations = kwargs.get('num_iterations', 1)
            log.info(f"Initialized with COUNT mode, num_iterations: {self.num_iterations}")
        # elif iteration_mode == IterationMode.HAMMING:
        #     self.hamming_distance = kwargs.get('hamming_distance', 5)
        #     log.info(f"Initialized with HAMMING mode, hamming_distance: {self.hamming_distance}")

        if aligner == Aligner.BWA:
            self.bwa_mem_args = kwargs.get("bwa_mem_args", "")

        # Make output path if it doesn't exist
        if not os.path.exists(output_path):
            os.makedirs(output_path)

    def run_sample(self, sample_id, read1_path, read2_path, ref_path: str):
        """
        Run the iterative alignment on a sample.
        """

        log.info("Running iterative alignment on sample: %s", sample_id)

        # Create execution directory in output folder
        execution_dir = os.path.join(self.output_path, sample_id)
        if not os.path.exists(execution_dir):
            os.makedirs(execution_dir)

        # Setup
        previous_ref_path = ref_path

        # Iterate over the number of iterations
        for i in range(1, self.max_iterations + 1):
            log.info(f"Iteration {i}")

            # Create output directory for this iteration in the execution directory
            iteration_dir = os.path.join(execution_dir, f"iteration_{i}")
            if not os.path.exists(iteration_dir):
                os.makedirs(iteration_dir)

            # Create index dir
            index_dir = os.path.join(iteration_dir, "index")
            if not os.path.exists(index_dir):
                os.makedirs(index_dir)

            # Copy the previous reference genome to the iteration directory
            current_ref_path = os.path.join(index_dir, f"ref_iter_{i}.fasta")
            shutil.copy(previous_ref_path, current_ref_path)

            # Align the sequences to the reference genome
            self.align(sample_id, i, iteration_dir, read1_path, read2_path, current_ref_path)

            # # Call variants
            # self.call_variants(i, sample_id, execution_dir)

            # # Mask the reference genome with the variants
            # self.mask(i, sample_id, execution_dir, ref_path)


    def align(self, sample_id: str, iter_num: int, iteration_dir: str, read1_path: str, read2_path: str, ref_path: str):
        """
        Align the sequences to the reference genome.
        """
        log.info("Aligning")

        # Switch on aligner
        if self.aligner == Aligner.BWA:
            # Call BWA index
            LogSubprocess().p_open(['bwa', 'index', ref_path])

            # Align
            log_subprocess = LogSubprocess()
            bwa_mem_cmd = ["bwa", "mem", "-t", str(self.num_cores), ref_path, read1_path, read2_path]
            bwa_mem_proc = log_subprocess.p_open(bwa_mem_cmd, stdout=subprocess.PIPE)

            # Convert to BAM
            samtools_view_to_bam_cmd = ["samtools", "view", "-@", str(self.num_cores), "-Sb", "-"]
            samtools_view_proc = log_subprocess.p_open(samtools_view_to_bam_cmd, stdin=bwa_mem_proc.stdout, stdout=subprocess.PIPE)
            bwa_mem_proc.stdout.close()

            # Sort
            samtools_sort_cmd = ["samtools", "sort", "-@", str(self.num_cores), "-o", os.path.join(iteration_dir, f"{sample_id}_iter_{iter_num}.bam"), "-"]
            samtools_sort_proc = log_subprocess.p_open(samtools_sort_cmd, stdin=samtools_view_proc.stdout)
            samtools_view_proc.stdout.close()

            # Error checking
            bwa_mem_proc.check_return_code()
            samtools_view_proc.check_return_code()
            samtools_sort_proc.check_return_code()


# per sample, per segment, per iteration

# 1. Align the sequences to the reference genome
# 2. Call variants
# 3. Mask the reference genome with the variants"""


# define enum for alignmer to use
# define enum for variant caller to use
# define enum for masking strategy