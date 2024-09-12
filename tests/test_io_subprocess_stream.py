"""
Subprocess stream tests
"""

# pylint: disable=missing-function-docstring,missing-class-docstring

import os
import unittest

from crick_genome_tools.io.subprocess_stream import SubprocessStream
from tests.utils import with_temporary_folder


CONTENT = "test_content"
READ_NAME = "@NB501505:171:H3KMGAFX3:1:21208:17616:17963 1:N:0:AGATCTCGGT"


class TestBuildFastaFromTopHits(unittest.TestCase):

    def test_subprocess_stream_gzip_read(self):
        """Test with read gzip file"""

        # Test and Assert
        with SubprocessStream(["gzip", "-c", "-d", "tests/data/io/fastq/sc_10k.fastq.gz"], mode="r") as stream:
            for idx, line in enumerate(stream):
                if idx == 0:
                    line_str = line.decode("UTF-8").strip()
                    self.assertEqual(line_str, READ_NAME)

    @with_temporary_folder
    def test_subprocess_stream_gzip_write(self, tmp_path):
        """Test with write gzip file"""

        # Setup
        filename = os.path.join(tmp_path, "test.gz")
        with open(filename, "w", encoding="UTF-8") as file:
            stream = SubprocessStream(["gzip", "-c"], stdout=file, mode="w")

            # Test
            stream.write(CONTENT.encode("UTF-8"))
            stream.close()

        # Assert
        stream = SubprocessStream(["gzip", "-c", "-d", filename], mode="r")
        assert stream is not None

        for idx, line in enumerate(stream):
            if idx == 0:
                self.assertEqual(line.decode("UTF-8"), CONTENT)
