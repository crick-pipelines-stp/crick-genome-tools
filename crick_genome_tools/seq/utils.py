"""
Utility functions for manipulating sequences.
"""


def hamming_distance(seq1, seq2):
    """
    Calculate the Hamming distance between two sequences.
    The Hamming distance is the number of positions at which the corresponding characters are different.
    If the sequences are of different lengths, the Hamming distance is calculated as the sum of:
    - The Hamming distance for the common length
    - The absolute difference in length between the two sequences

    Args:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.

    Returns:
        int: The Hamming distance between the two sequences.
    """
    # Calculate the base Hamming distance for the common length
    min_len = min(len(seq1), len(seq2))
    distance = sum(c1 != c2 for c1, c2 in zip(seq1[:min_len], seq2[:min_len]))

    # Add the penalty for length difference
    distance += abs(len(seq1) - len(seq2))

    return distance
