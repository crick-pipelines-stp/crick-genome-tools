# pylint: disable=missing-function-docstring,missing-class-docstring


import os

import pytest
from assertpy import assert_that

from crick_genome_tools.contig import filter_viral_contigs


class TestContigs:
    def test_custom(self):
        filter_viral_contigs.filter_viral_contigs(50, 0.5, "PB2,PB1,PA,HA,NP,NA,MP,NS", "test_data.tsv")
        raise Exception("finish me")