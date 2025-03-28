"""
Tests for reading and writing fastq files.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring

from assertpy import assert_that

from crick_genome_tools.io.fastq_file import FastqFile

class TestFastqFile:
    def test_io_fastq_file_read_uncompressed(self):
        # Setup
        fastq_file = "tests/data/io/fastq/test.fastq"
        fastq = FastqFile(fastq_file)

        # Test and assert
        for name, seq, qual in fastq.open_read_iterator(as_string=True):
            assert_that(name).is_equal_to("NB501505:171:H3KMGAFX3:4:21612:13641:20150 1:N:0:CTATAGTCTT")
            assert_that(seq).is_equal_to("AGCAGATAGTGAGGAAAGTTGAGCCAATAATGACGTGAAGTCCGTGGAAGCC")
            assert_that(qual).is_equal_to("AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEA")
            break

    def test_io_fastq_file_read_compressed(self):
        # Setup
        fastq_file = "tests/data/io/fastq/sc_10k.fastq.gz"
        fastq = FastqFile(fastq_file)

        # Test and assert
        idx = 0
        for name, _, _ in fastq.open_read_iterator(as_string=True):
            if idx == 0:
                assert_that(name).is_equal_to("NB501505:171:H3KMGAFX3:1:21208:17616:17963 1:N:0:AGATCTCGGT")
            idx += 1  # Avoid early closure errors

    def test_io_fastq_write_file(self, tmp_path):
        # Setup
        fastq_file = "tests/data/io/fastq/test.fastq"
        fastq_read = FastqFile(fastq_file)
        fastq_write = FastqFile(str(tmp_path / "test.fastq.gz"))

        # Test
        write_stream = fastq_write.open_write_stream()
        for name, seq, qual in fastq_read.open_read_iterator(as_string=True):
            FastqFile.write_read(write_stream, name, seq, qual)
        write_stream.close()

        # Assert
        output_file = str(tmp_path / "test.fastq.gz")
        fastq = FastqFile(output_file)
        for name, seq, qual in fastq.open_read_iterator(as_string=True):
            assert_that(name).is_equal_to("NB501505:171:H3KMGAFX3:4:21612:13641:20150 1:N:0:CTATAGTCTT")
            assert_that(seq).is_equal_to("AGCAGATAGTGAGGAAAGTTGAGCCAATAATGACGTGAAGTCCGTGGAAGCC")
            assert_that(qual).is_equal_to("AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEA")
            break
