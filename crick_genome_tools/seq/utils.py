"""
Utility functions for manipulating sequences.
"""

def rev_comp(seq):
    """Reverse complement a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement[base] for base in reversed(seq))


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


def cumulative_hamming_distance(dict1, dict2):
    """
    Calculate the cumulative Hamming distance between two dictionaries of sequences.

    The function compares sequences stored in two dictionaries (dict1 and dict2) by their keys.
    If a key exists in only one dictionary, the sequence from the other dictionary is considered
    an empty string, and the total distance for that entry is the full length of the sequence.
    The Hamming distance is calculated for each matching entry and summed up to get the total.

    Args:
        dict1 (dict): A dictionary where keys are sequence identifiers and values are sequences (strings).
        dict2 (dict): A second dictionary with the same structure as dict1.

    Returns:
        int: The total cumulative Hamming distance between all matching entries in the two dictionaries.
    """

    total_distance = 0
    all_keys = set(dict1.keys()).union(set(dict2.keys()))  # Collect all unique keys

    for key in all_keys:
        seq1 = dict1.get(key, "")  # Get sequence from dict1 or empty string if missing
        seq2 = dict2.get(key, "")  # Get sequence from dict2 or empty string if missing
        total_distance += hamming_distance(seq1, seq2)

    return total_distance
