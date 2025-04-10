"""
Test the SubprocessStream class.
"""

# pylint: disable=missing-function-docstring,missing-class-docstring

from assertpy import assert_that

from crick_genome_tools.io.subprocess_stream import SubprocessStream


CONTENT = "test_content"
READ_NAME = "@NB501505:171:H3KMGAFX3:1:21208:17616:17963 1:N:0:AGATCTCGGT"


class TestSubprocessStream:
    def test_io_subproc_gzip_read(self):
        with SubprocessStream(["gzip", "-c", "-d", "tests/data/io/fastq/sc_10k.fastq.gz"], mode="r") as stream:
            for idx, line in enumerate(stream):
                if idx == 0:
                    line_str = line.decode("utf-8").strip()
                    assert_that(line_str).is_equal_to(READ_NAME)

    def test_io_subproc_gzip_write_and_read_back(self, tmp_path):
        out_file = tmp_path / "test.gz"

        # Write
        with open(out_file, "w", encoding="utf-8") as f:
            stream = SubprocessStream(["gzip", "-c"], stdout=f, mode="w")
            stream.write(CONTENT.encode("utf-8"))
            stream.close()

        # Read
        with SubprocessStream(["gzip", "-c", "-d", str(out_file)], mode="r") as stream:
            for idx, line in enumerate(stream):
                if idx == 0:
                    assert_that(line.decode("utf-8")).is_equal_to(CONTENT)
