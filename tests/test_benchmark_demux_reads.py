"""
Benchmarking tests for generate reads
"""

import io
import os
import random
import tempfile

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
    # fastq_file = os.path.join("tests/data/seq/L002_R1.fastq")
    # sample_dict = {
    #     "sample_1": "GCTT,NTAT",
    #     "sample_2": "ACGT,AGGT",
    # }

    fastq_file = os.path.join("tests/data/seq/sub_read_L002_R1.fastq")
    # sample_dict = {
    #     "sample_1": "ACTGGTGTCG-CAAGTCCTGT",
    #     "sample_2": "GCAGTCTTAT-CCGGCCATTA",
    #     "sample_3": "CGCAGAACTT-GACGTCGATA",
    #     "sample_4": "GCCTAGGACT-AGCCGTTCTC",
    #     "sample_5": "ATGGTTGACT-TGGCACAAGC",
    #     "sample_6": "GTGTTATCTC-AGTCTGGTGT",
    # }
    sample_dict = {
        "sample_1": "ACTGGTGTCG-CAAGTCCTGT",
        "sample_2": "AGGTGGCTAC+CCACGTAACG",
        "sample_3": "TATCACTCTC+ACCTTGTTCT",
        "sample_4": "AGGTGGCTAC+CCACGTAACG",
        "sample_5": "ATGGTTGACT-TGGCACAAGC",
        "sample_6": "GCAGTCTTAT+GGGGGGGGGG",
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


def generate_sample_dict(num_samples=10):
    """
    Generate a dictionary of sample names mapped to random dual-index barcodes.

    Args:
        num_samples (int): Number of samples to generate.

    Returns:
        dict: A dictionary where keys are sample names and values are barcodes.
    """
    sample_dict = {}
    for i in range(1, num_samples + 1):
        # Use + or - to simulate dual-index formats
        separator = "+"
        barcode1 = "".join(random.choices("ACGT", k=10))
        barcode2 = "".join(random.choices("ACGT", k=10))
        sample_dict[f"sample_{i}"] = f"{barcode1}{separator}{barcode2}"
    return sample_dict


def generate_barcode_fastq(sample_dict, num_reads=1_000_000, read_length=100):
    """
    Generate a FASTQ file-like object with reads containing barcodes from the sample dictionary.

    Args:
        sample_dict (dict): Dictionary of sample names and barcodes.
        num_reads (int): Number of reads to generate.
        read_length (int): Length of each read.

    Returns:
        io.StringIO: In-memory file-like object containing the FASTQ content.
    """
    barcodes = list(sample_dict.values())
    bases = ["A", "C", "G", "T"]
    qualities = ["I"] * read_length
    fastq_lines = []

    for i in range(num_reads):
        barcode = random.choice(barcodes)
        header = f"@LH00442{i} 1:N:0:{barcode}"
        sequence = "".join(random.choices(bases, k=read_length))
        plus_line = "+"
        quality = "".join(qualities)
        fastq_lines.extend([header, sequence, plus_line, quality])

    fastq_content = "\n".join(fastq_lines)
    return io.StringIO(fastq_content)  # In-memory file-like object


@pytest.mark.benchmark(group="demux-reads-illumina", min_rounds=5)
@pytest.mark.only_run_with_direct_target
def test_benchmarking_1M_reads_demux_fastq(*args, tmp_path, benchmark):  # pylint: disable=unused-argument
    """Test generate reads to a fastq file"""

    # Variables
    max_time_ms = int(os.getenv("MAX_BENCHMARK_GENREADS_MS", "200"))

    # Setup
    # Generate the FASTQ content in memory
    # fastq_content = generate_random_fastq()
    # fastq_file = io.StringIO(fastq_content)

    sample_dict = generate_sample_dict(num_samples=10)

    max_hamming_distance = 1
    output_path = tmp_path


    # Generate FASTQ content
    fastq_content = generate_barcode_fastq(sample_dict, num_reads=1_000_000).getvalue()

    # Write to a temporary file
    with tempfile.NamedTemporaryFile(mode='w+', delete=True, suffix='.fastq') as temp_file:
        temp_file.write(fastq_content)
        temp_file.flush()  # Ensure data is written to disk

        # Now pass the file path
        # demultiplex_fastq_by_barcode(temp_file.name, sample_dict, max_hamming_distance, output_path)
        benchmark(lambda: (demultiplex_fastq_by_barcode(temp_file.name, sample_dict, max_hamming_distance, output_path)))


    # Test
    # benchmark(lambda: (demultiplex_fastq_by_barcode(fastq_file, sample_dict, max_hamming_distance, output_path)))

    # Assert
        mean_ms = benchmark.stats["mean"] * 1000
        if mean_ms > max_time_ms:
            pytest.fail(f"Benchmark failed: expected - {max_time_ms}ms | actual - {round(mean_ms, 2)}ms")

