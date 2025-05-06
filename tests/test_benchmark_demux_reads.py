"""
Benchmarking tests for generate reads
"""

import os

import pytest

# from .test_seq_barcode_demux import test_demultiplex_fastq_by_barcode_valid
from crick_genome_tools.seq.barcode_demux import demultiplex_fastq_by_barcode


# TEST_CONFIG_SOURCE_FOLDER = "tests/data/config"
# TEST_FASTA_SOURCE_FOLDER = "tests/data/refs"


@pytest.mark.benchmark(group="demux-reads-illumina", min_rounds=5)
@pytest.mark.only_run_with_direct_target
def test_benchmarking_generate_reads_fastq(*args, tmp_path, benchmark):  # pylint: disable=unused-argument
    """Test generate reads to a fastq file"""

    # Variables
    max_time_ms = int(os.getenv("MAX_BENCHMARK_GENREADS_MS", "200"))

    # Setup
    fastq_file = os.path.join("tests/data/seq/L002_R1.fastq")
    sample_dict = {
        "sample_1": "GCTT,NTAT",
        "sample_2": "ACGT,AGGT",
    }
    max_hamming_distance = 1
    output_path = tmp_path

    # Test
    benchmark(lambda: (demultiplex_fastq_by_barcode(fastq_file, sample_dict, max_hamming_distance, output_path)))

    # Assert
    mean_ms = benchmark.stats["mean"] * 1000
    if mean_ms > max_time_ms:
        pytest.fail(f"Benchmark failed: expected - {max_time_ms}ms | actual - {round(mean_ms, 2)}ms")

    # raise ValueError
