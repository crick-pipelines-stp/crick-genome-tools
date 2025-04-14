"""
Benchmarking tests for generate reads
"""

import os

import pytest

from .test_seq_barcode_demux import gen_index_and_hamming_dict


# TEST_CONFIG_SOURCE_FOLDER = "tests/data/config"
# TEST_FASTA_SOURCE_FOLDER = "tests/data/refs"


@pytest.mark.benchmark(group="generate-reads-illumina", min_rounds=5)
@with_temporary_folder
@pytest.mark.only_run_with_direct_target
def test_benchmarking_generate_reads_fastq_pe(*args, tmp_path, benchmark):  # pylint: disable=unused-argument
    """Test generate reads to a fastq file"""

    # Variables
    num_reads = 10000
    max_time_ms = int(os.getenv("MAX_BENCHMARK_GENREADS_MS", "200"))

    # Setup
    profile = "illumina_pe_nextera"
    # fasta_dir = os.path.join(TEST_FASTA_SOURCE_FOLDER, "valid")
    # config_path = os.path.join(TEST_CONFIG_SOURCE_FOLDER, "full_config.yml")
    # output_path = os.path.join(tmp_path, "output.fastq.gz")
    # test = GenerateReads(config_path, fasta_dir)
    fastq_file = os.path.join("tests/data/io/fastq/test.fastq")
    max_hamming_distance = 1
    output_path = tmp_path
    test = gen_index_and_hamming_dict(fastq_file, dict, max_hamming_distance)

    # Test
    # benchmark(test.run, profile, output_path, num_reads)
    benchmark(test.run, profile, output_path)

    # # Assert
    # mean_ms = benchmark.stats["mean"] * 1000
    # if mean_ms > max_time_ms:
    #     pytest.fail(f"Benchmark failed: expected - {max_time_ms}ms | actual - {round(mean_ms, 2)}ms")