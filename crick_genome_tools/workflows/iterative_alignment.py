"""
Performs iterative alignment of sequences to a reference genome.
"""

import csv
import logging
import os
import shutil
from enum import Enum

from crick_genome_tools.io.command_chain import CommandChain
from crick_genome_tools.io.fasta import Fasta
from crick_genome_tools.io.log_subprocess import LogSubprocess
from crick_genome_tools.io.samtools import parse_flagstat
from crick_genome_tools.seq.utils import cumulative_hamming_distance


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
        # if iteration_mode == IterationMode.COUNT:
        #     self.num_iterations = kwargs.get("num_iterations", 1)
        #     log.info(f"Initialized with COUNT mode, num_iterations: {self.num_iterations}")
        # elif iteration_mode == IterationMode.HAMMING:
        #     self.hamming_distance = kwargs.get('hamming_distance', 5)
        #     log.info(f"Initialized with HAMMING mode, hamming_distance: {self.hamming_distance}")

        # Aligner specific configurations
        if aligner == Aligner.BWA:
            # Initialize dynamic params for BWA (mem, mmpen, gappen)
            self.aligner_params = {
                "bwa_args": kwargs.get("bwa_args", []),
                "mem": kwargs.get("mem_start", 18),
                "mmpen": kwargs.get("mmpen_start", 10),
                "gappen": kwargs.get("gappen_start", 5),
            }
            self.aligner_param_increments = {
                "mem": kwargs.get("mem_increment", 2),
                "mmpen": kwargs.get("mmpen_increment", 1),
                "gappen": kwargs.get("gappen_increment", 1),
            }
            self.aligner_param_endpoints = {
                "mem": kwargs.get("mem_end", 30),
                "mmpen": kwargs.get("mmpen_end", 15),
                "gappen": kwargs.get("gappen_end", 10),
            }

            # Log these params
            log.info(
                f"Initialized with BWA aligner, params: {self.aligner_params}, increments: {self.aligner_param_increments}, endpoints: {self.aligner_param_endpoints}"
            )

        # Init other params
        # self.num_iterations = kwargs.get("realignment_minqscore", 1)

    def run_sample(self, sample_id: str, read1_path: str, read2_path: str, ref_path: str):
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
        iteration_metrics = []

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

            # Create logs dir
            log_dir = os.path.join(iteration_dir, "logs")
            if not os.path.exists(log_dir):
                os.makedirs(log_dir)

            # Copy the previous reference genome to the iteration directory
            current_ref_path = os.path.join(index_dir, f"ref_iter_{i}.fasta")
            shutil.copy(previous_ref_path, current_ref_path)

            # Align the sequences to the reference genome
            bam_file, flagstat_file = self.align(sample_id, i, iteration_dir, read1_path, read2_path, current_ref_path, log_dir)

            # Call variants
            realigned_bam_file, realigned_bai_file = self.realign_sequences(sample_id, i, iteration_dir, current_ref_path, log_dir, bam_file)

            # Generate a consensus sequence
            consensus_fasta_path = self.gen_consesnsus(sample_id, i, iteration_dir, current_ref_path, log_dir, realigned_bam_file)

            # Parse the flagstat file
            total_reads, mapped_reads, alignment_rate = parse_flagstat(flagstat_file)

            # Check diff
            prev_fasta = Fasta.read_fasta_file(previous_ref_path)
            con_fasta = Fasta.read_fasta_file(consensus_fasta_path)
            hamming_dist = cumulative_hamming_distance(prev_fasta, con_fasta)
            log.info(f"Hamming distance: {hamming_dist}")

            # Update the alignment parameters
            log.info("Updating params")
            self.update_params(i)

            # Store metrics for this iteration
            metrics = {
                "Sample Id": sample_id,
                "Iteration": i,
                "Total Reads": total_reads,
                "Mapped Reads": mapped_reads,
                "Alignment Rate (%)": alignment_rate,
                "Hamming Distance": hamming_dist,
            }
            iteration_metrics.append(metrics)

            # Update the reference path for the next iteration
            previous_ref_path = consensus_fasta_path

        log.info("Finished iterative alignment")

        # Copy the final consensus sequence to the output directory
        final_consensus_path = os.path.join(execution_dir, f"{sample_id}_final_consensus.fasta")
        shutil.copy(consensus_fasta_path, final_consensus_path)

        # Copy the final bam and bai files to the output directory
        final_bam_path = os.path.join(execution_dir, f"{sample_id}_final.bam")
        final_bai_path = os.path.join(execution_dir, f"{sample_id}_final.bai")
        shutil.copy(realigned_bam_file, final_bam_path)
        shutil.copy(realigned_bai_file, final_bai_path)

        # Copy the final flagstat file to the output directory
        final_flagstat_path = os.path.join(execution_dir, f"{sample_id}_final.flagstat")
        shutil.copy(flagstat_file, final_flagstat_path)

        # Output the iteration metrics to a CSV file
        metrics_file = os.path.join(execution_dir, f"{sample_id}_iteration_metrics.csv")
        with open(metrics_file, "w", newline="", encoding="UTF-8") as csvfile:
            fieldnames = ["Sample Id", "Iteration", "Total Reads", "Mapped Reads", "Alignment Rate (%)", "Hamming Distance"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()
            for data in iteration_metrics:
                writer.writerow(data)

    def align(self, sample_id: str, iter_num: int, iteration_dir: str, read1_path: str, read2_path: str, ref_path: str, log_dir: str):
        """
        Align the sequences to the reference genome.
        """
        log.info("Aligning")

        # Init file names
        bam_file = os.path.join(iteration_dir, f"{sample_id}_iter_{iter_num}.bam")
        flagstat_file = os.path.join(iteration_dir, f"{sample_id}_iter_{iter_num}.flagstat")

        # Switch on aligner
        if self.aligner == Aligner.BWA:
            # Call BWA index
            CommandChain.command_to_logfile(["bwa-mem2", "index", ref_path], os.path.join(log_dir, f"{sample_id}_iter_{iter_num}.refindex.log"))

            # Define the BWA mem command using dynamic params
            bwa_command = (
                [
                    "bwa-mem2",
                    "mem",
                    "-t",
                    str(self.num_cores),
                    "-R",
                    f"@RG\tID:{sample_id}\tSM:{sample_id}\tLB:{sample_id}\tPL:ILLUMINA",
                    "-k",
                    str(self.aligner_params["mem"]),
                    "-B",
                    str(self.aligner_params["mmpen"]),
                    "-O",
                    str(self.aligner_params["gappen"]),
                ]
                + self.aligner_params["bwa_args"]
                + [ref_path, read1_path, read2_path]
            )
            log.info(f"Running BWA mem with command: {bwa_command}")

            # Define the alignment command chain and run
            commands = [
                bwa_command,  # BWA mem
                ["samtools", "view", "-@", str(self.num_cores), "-Sb", "-"],  # Convert to BAM
                ["samtools", "sort", "-@", str(self.num_cores), "-o", bam_file, "-"],  # Sort
            ]
            command_chain = CommandChain(commands)
            command_chain.run()

        # Index the bam file
        LogSubprocess().p_open(["samtools", "index", "-@", str(self.num_cores), bam_file]).check_return_code()

        # Run samtools flagstat
        CommandChain.command_to_file(["samtools", "flagstat", "-@", str(self.num_cores), bam_file], flagstat_file)

        return bam_file, flagstat_file

    def realign_sequences(self, sample_id: str, iter_num: int, iteration_dir: str, ref_path: str, log_dir: str, bam_file: str):
        """
        Call variants on the aligned sequences.
        """
        log.info("Realigning")

        # Init file names
        realigned_bam_file = os.path.join(iteration_dir, f"{sample_id}_iter_{iter_num}.realigned.bam")
        realigned_bai_file = os.path.join(iteration_dir, f"{sample_id}_iter_{iter_num}.realigned.bai")

        # Realign bam with abra2
        CommandChain.command_to_logfile(
            ["abra2", "--threads", str(self.num_cores), "--in", bam_file, "--out", realigned_bam_file, "--ref", ref_path, "--index"],
            os.path.join(log_dir, f"{sample_id}_iter_{iter_num}.realign.log"),
        )

        return realigned_bam_file, realigned_bai_file

    def gen_consesnsus(self, sample_id: str, iter_num: int, iteration_dir: str, ref_path: str, log_dir: str, bam_file: str):
        """
        Generate a consensus sequence from the aligned sequences.
        """
        log.info("Generating consensus")

        # Init file names
        vcf_file = os.path.join(iteration_dir, f"{sample_id}_iter_{iter_num}.vcf.gz")
        fasta_file = os.path.join(iteration_dir, f"{sample_id}_iter_{iter_num}.consensus.fasta")

        # Call variants
        commands = [
            ["bcftools", "mpileup", "--threads", str(self.num_cores), "-Q 20", "-L 10000", "-Ep", "-f", ref_path, bam_file],  # pileup
            ["bcftools", "call", "--threads", str(self.num_cores), "-c", "-Oz", "-"],  # call
        ]
        command_chain = CommandChain(commands, output_file=vcf_file)
        command_chain.run()

        # Index the vcf file
        LogSubprocess().p_open(["bcftools", "index", vcf_file]).check_return_code()

        # Generate consensus
        CommandChain.command_to_logfile(
            ["bcftools", "consensus", "-f", ref_path, "-o", fasta_file, vcf_file], os.path.join(log_dir, f"{sample_id}_iter_{iter_num}.consensus.log")
        )

        return fasta_file

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
            if endpoint is not None and (
                (increment > 0 and self.aligner_params[param] > endpoint) or (increment < 0 and self.aligner_params[param] < endpoint)
            ):
                self.aligner_params[param] = endpoint
