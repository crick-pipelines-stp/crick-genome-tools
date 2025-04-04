
import pytest
from assertpy import assert_that


from crick_genome_tools.seq.barcode_demux import test_function

class TestBarcodeDemux:
    def test_barcode_demux_function(self):
        assert_that(test_function()).is_true()