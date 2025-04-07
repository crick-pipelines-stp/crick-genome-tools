# """
# Tests itr mask
# """

# # pylint: disable=missing-function-docstring,missing-class-docstring

# import os
# from assertpy import assert_that
# from crick_genome_tools.seq.mask_itr import get_variable_itr_regions

# from crick_genome_tools.io.fasta import Fasta

# FIVE_PRIME_FLIP_ITR = "GCGCGCGCCCGGGCAAAGCCCGGGCGTCGGGCGACCTTTGGTCGCCCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCA"
# THREE_PRIME_FLIP_ITR = "GCGCGCGCTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGG"

# class TestMaskItr:
#     def test_itr(self):
#         # Setup
#         fasta = Fasta.read_fasta_file("tests/data/seq/aav.fasta")
#         itr_seq1 = fasta["P220_AAV_2100bp"][0:200]
#         result = get_variable_itr_regions(itr_seq1)
#         print(result)
#         raise NotImplementedError
