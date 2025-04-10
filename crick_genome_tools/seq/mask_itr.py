"""
Mask the variable ITR regions.

AAV transgene cassettes are flanked by ITR regions, which can exist in flip or flop
configurations determined by the order of B'/B and C'/C palindromic elements.
As there are two ITR regions in the transgene cassette, there are a total of 4 states
that can naturally occur.

In 'flip' orientation, the C'-C palindromes are located nearest to the 5' end of the ITR
In 'flop' orientation, to B'-B palindromes are located nearest to the 5' end of the ITR

In this code, B' and B elements are referred to as b_prime and b (similarly for C)
"""

import logging

from crick_genome_tools.seq.utils import rev_comp


log = logging.getLogger(__name__)


def get_variable_itr_regions(itr_seq):
    """
    Identify ITR regions.

    :param str itr_seq: a single ITR sequence
    :returns:
        - locs_in_itr: the varaible positions within the ITR


        The two ITR configurations are shown below for the subregions we are
        interested in.

            flip orientation
            5' ... C' AAA_loop - C - T(single nt) - B' - TTT_loop - B ... 3'

            flop orientation
            5' ... B' AAA_loop - B - A(single nt) - C' - TTT_loop - C ... 3'

        The B-B' and C-C' elements form a T-shaped hairpin structure.
        Each end of the hairpin contains a non-palindromic 3 nt loop. The location of
        the B-C region is identified by looking for the two loop structures,
        which should be 19 nt apart.
        The AAA lop should be nearest to the 5' end and the TTT loop nearest to the
        3' end.
    """
    # The location of the B-C palindromes is calculated in relation to the
    # AAA and TTT hairpin loops. These must exist, and be separated by 19 nt.
    aaa_loop_index = itr_seq.find("AAA")  # Always closer to 5'
    ttt_loop_index = itr_seq.find("TTT")  # Always closer to 3'
    loop_distance = ttt_loop_index - aaa_loop_index - 3
    log.info(f"AAA loop found at position {aaa_loop_index}")
    log.info(f"TTT loop found at position {ttt_loop_index}")
    log.info(f"Loop distance is {loop_distance} nt")

    # Do some checks to ensure correct transgene plasmid used.
    if -1 in [aaa_loop_index, ttt_loop_index]:
        raise ValueError(
            "ITR sequence missing 'AAA or TTT' hairpin sequence." "Please check your transgene plasmid sequence and ITR location parameters"
        )
    if ttt_loop_index - aaa_loop_index - 3 != 19:
        raise ValueError(
            "Distance between 'TTT' loop and AAA should be 19 nt " "Please check your transgene plasmid sequence and ITR location parameters"
        )

    # Extract all the B and C elements in order they are in the plasmid.
    elem1 = itr_seq[aaa_loop_index - 9 : aaa_loop_index]
    elem2 = itr_seq[aaa_loop_index + 3 : aaa_loop_index + 12]
    elem3 = itr_seq[ttt_loop_index - 9 : ttt_loop_index]
    elem4 = itr_seq[ttt_loop_index + 3 : ttt_loop_index + 12]
    base_between_b_c = itr_seq[aaa_loop_index + 12]

    if base_between_b_c == "T":  # Flip
        e1 = "C'"
        e2 = "C"
    elif base_between_b_c == "A":  # Flop
        e1 = "B'"
        e2 = "B"
    else:
        raise ValueError(
            f"Base between B and C elements should be 'T' or 'A' not {base_between_b_c}"
            " Please check your transgene plasmid sequence and ITR location parameters"
        )

    # Check if the C',C and B',B elements are palindromes. Warn if not
    if elem1 != rev_comp(elem2):
        raise ValueError(f"{e1} and {e2} are not palindromic - " "Please check your transgene plasmid sequence and ITR location parameters")
    if elem3 != rev_comp(elem4):
        raise ValueError(f"{e1} and {e2} are not palindromic - " "Please check your transgene plasmid sequence and ITR location parameters")

    # We don't need to know whether this ITR is in flip or flop orientation
    # Just need to reorder the elements (B', B, C' C), to get both.
    itr_plasmid_orientation = f"{elem1}AAA{elem2}T{elem3}TTT{elem4}"
    itr_alternative_orientation = f"{elem3}AAA{elem4}A{elem1}TTT{elem2}"

    # Now find variable locations
    elem1_pos = aaa_loop_index - 9
    locs_in_itr = [
        i + elem1_pos
        for i, (original, alternative) in enumerate(zip(itr_plasmid_orientation, itr_alternative_orientation))
        if original != alternative
    ]

    return locs_in_itr
