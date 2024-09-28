"""
Performs iterative alignment of sequences to a reference genome.
"""

import logging
import os
import shutil
from enum import Enum

from crick_genome_tools.io.log_subprocess import LogSubprocess
from crick_genome_tools.io.command_chain import CommandChain


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

    def __init__(
        self, num_cores: int, output_path: str, min_iterations: int, max_iterations: int, aligner: Aligner, iteration_mode: IterationMode, **kwargs
    ):
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
            self.num_iterations = kwargs.get("num_iterations", 1)
            log.info(f"Initialized with COUNT mode, num_iterations: {self.num_iterations}")
        # elif iteration_mode == IterationMode.HAMMING:
        #     self.hamming_distance = kwargs.get('hamming_distance', 5)
        #     log.info(f"Initialized with HAMMING mode, hamming_distance: {self.hamming_distance}")

        # Aligner specific configurations
        if aligner == Aligner.BWA:
            # Initialize dynamic params for BWA (mem, mmpen, gappen)
            self.aligner_params = {
                'bwa_args': kwargs.get('bwa_args', []),
                'mem': kwargs.get('mem_start', 20),
                'mmpen': kwargs.get('mmpen_start', 10),
                'gappen': kwargs.get('gappen_start', 5)
            }
            self.aligner_param_increments = {
                'mem': kwargs.get('mem_increment', 0),
                'mmpen': kwargs.get('mmpen_increment', 0),
                'gappen': kwargs.get('gappen_decrement', 0)
            }
            self.aligner_param_endpoints = {
                'mem': kwargs.get('mem_end', None),
                'mmpen': kwargs.get('mmpen_end', None),
                'gappen': kwargs.get('gappen_end', None)
            }

            # Log these params
            log.info(f"Initialized with BWA aligner, params: {self.aligner_params}, increments: {self.aligner_param_increments}, endpoints: {self.aligner_param_endpoints}")

    def run_sample(self, sample_id, read1_path, read2_path, ref_path: str):
        """
        Run the iterative alignment on a sample.
        """

        log.info("Running iterative alignment on sample: %s", sample_id)

        # Make output path if it doesn't exist
        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path)

        # Create execution directory in output folder
        execution_dir = os.path.join(self.output_path, sample_id)
        if not os.path.exists(execution_dir):
            os.makedirs(execution_dir)

        # Setup
        previous_ref_path = ref_path

        # Iterate over the number of iterations
        for i in range(1, self.max_iterations + 1):
            log.info(f"Iteration {i}")

            # Create output directory for this iteration in the execution directory
            iteration_dir = os.path.join(execution_dir, f"iteration_{i}")
            if not os.path.exists(iteration_dir):
                os.makedirs(iteration_dir)

            # Create index dir
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

        # Init file names
        bam_file = os.path.join(iteration_dir, f"{sample_id}_iter_{iter_num}.bam")

        # Switch on aligner
        if self.aligner == Aligner.BWA:
            # Call BWA index
            LogSubprocess().p_open(["bwa-mem2", "index", ref_path]).check_return_code()

            # Define the BWA mem command using dynamic params
            bwa_command = [
                "bwa-mem2", "mem",
                "-t", str(self.num_cores),
                "-R", f"@RG\tID:{sample_id}\tSM:{sample_id}\tLB:{sample_id}\tPL:ILLUMINA",
                "-k", str(self.aligner_params['mem']),
                "-B", str(self.aligner_params['mmpen']),
                "-O", str(self.aligner_params['gappen'])] + self.aligner_params['bwa_args'] + [
                ref_path, read1_path, read2_path
            ]
            log.info(f"Running BWA mem with command: {bwa_command}")

            # Define the alignment command chain and run
            commands = [
                bwa_command,  # BWA mem
                ["samtools", "view", "-@", str(self.num_cores), "-Sb", "-"],  # Convert to BAM
                ["samtools", "sort", "-@", str(self.num_cores), "-o", bam_file, "-"]  # Sort
            ]
            command_chain = CommandChain(commands)
            command_chain.run()

        # Index the bam file
        LogSubprocess().p_open(["samtools", "index", "-@", str(self.num_cores), bam_file]).check_return_code()

        # Run samtools flagstat
        CommandChain.command_to_file(["samtools", "flagstat", "-@", str(self.num_cores), bam_file], os.path.join(iteration_dir, f"{sample_id}_iter_{iter_num}.flagstat"))

    def update_params(self, iteration: int):
        """
        Update the alignment parameters based on the param_mode and iteration number.
        """
        for param, increment in self.aligner_param_increments.items():
            # Update the parameter value if increment is defined
            if increment != 0:
                self.aligner_params[param] += increment

            # Cap the parameter if an endpoint is defined
            endpoint = self.aligner_param_endpoints.get(param)
            if endpoint is not None and ((increment > 0 and self.aligner_params[param] > endpoint) or
                                            (increment < 0 and self.aligner_params[param] < endpoint)):
                self.aligner_params[param] = endpoint
