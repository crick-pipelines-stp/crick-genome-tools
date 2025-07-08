import gzip
import os
import re
from collections import defaultdict

import numpy as np
from pybktree import BKTree

from crick_genome_tools.io.fastq_file import FastqFile


def group_samples_by_index_length(sample_index_dict: dict) -> list:
    """
    Groups samples by the lengths of their barcode index components.

    This function processes a dictionary mapping sample names to barcode data,
    which may be a string or a dictionary with keys like "barcode", "index", or "index2".
    It parses these barcodes, splits them into components (if needed), and groups samples
    into categories based on the lengths of their first and optional second index component.

    Args:
        sample_index_dict (dict): A dictionary where keys are sample names and values are
            either strings representing barcodes or dictionaries containing barcode/index fields.

    Returns:
        dict: A dictionary where keys are tuples representing the lengths of index components
            (e.g., (4, 4)), and values are dictionaries mapping sample names to lists of index parts.

    Raises:
        TypeError: If the input is not a dictionary or if a sample value is of an unsupported type.
        ValueError: If a sample's barcode splits into more than two components.
    """
    # Check if the input is a dictionary
    if not isinstance(sample_index_dict, dict):
        raise TypeError(f"{sample_index_dict} must be a dictionary.")

    result = defaultdict(dict)
    # Calculate barcode lengths and organize samples by these lengths within `result`
    for sample, value in sample_index_dict.items():
        index = []

        if isinstance(value, str):
            # if sample has a string value, treat it as a barcode
            index = re.split(r"[^A-Za-z]+", value)

        elif isinstance(value, dict):
            # if sample has a dict value, look for "barcode" or "index" keys
            if "barcode" in value:
                index = re.split(r"[^A-Za-z]+", value["barcode"])
            elif "index" in value:
                if "index2" in value:
                    # if both index and index2 are present, use both entries as the barcode
                    index = [value["index"], value["index2"]]
                else:
                    # if only index is present, use it as the barcode
                    index = re.split(r"[^A-Za-z]+", value["index"])
        else:
            raise TypeError(f"Sample '{sample}' has an unsupported value type: {type(value).__name__}")

        # Ensure that there are 2 indexes maximum
        if len(index) > 2:
            raise ValueError(f"Sample '{sample}' has more than two index components: {index}")

        # Determine the lengths of the first and second index
        first_len = len(index[0])
        last_len = 0
        if len(index) == 2:
            last_len = len(index[1])

        key = (first_len, last_len)
        result[key][sample] = index

    return dict(result)


def hamming_distance(seq1, seq2) -> int:
    """Computes the Hamming distance between two sequences."""
    if seq1 is None or seq2 is None:
        raise ValueError("Input sequences cannot be None.")

    return np.count_nonzero(np.frombuffer(seq1.encode(), dtype="S1") != np.frombuffer(seq2.encode(), dtype="S1"))


def crosscheck_barcode_proximity(barcodes: dict) -> list:
    """
    This function performs an all-vs-all comparison of barcode strings provided
    in the input dictionary. For each pair of barcodes, it calculates the Hamming
    distance and includes the pair in the output if the distance is within the
    specified threshold.

    Args:
        barcodes (dict): A dictionary mapping sample names to barcode strings.

    Returns:
        list: A list of tuples, where each tuple contains:
            - The first barcode string.
            - The second barcode string.
            - The Hamming distance between the two barcodes.

    Raises:
        ValueError: If any pair of barcodes has unequal lengths, as Hamming distance
                    requires strings of the same length.
    """
    if not isinstance(barcodes, dict):
        raise ValueError(f"{barcodes} must be a dictionary.")

    items = list(barcodes.items())
    similar_pairs = []

    for i in range(len(items)):
        sample1, bc1 = items[i]
        for j in range(i + 1, len(items)):
            sample2, bc2 = items[j]

            if len(bc1) != len(bc2):
                raise ValueError(f"Barcodes '{sample1}' and '{sample2}' are of unequal length")

            # Fast inline hamming distance with early exit
            dist = 0
            for a, b in zip(bc1, bc2):
                if a != b:
                    dist += 1
            similar_pairs.append((bc1, bc2, dist))
    return similar_pairs


def find_min_hamming_distances(grouped_sample_hamming_by_length: dict) -> dict:
    """
    Finds the minimum Hamming distance for each group length.

    Args:
        grouped_sample_hamming_by_length (dict): A dictionary where each key is a group length
            and the value is a set of tuples (barcode1, barcode2, hamming_distance).

    Returns:
        dict: A dictionary mapping each group length to the smallest Hamming distance found.
    """
    if not isinstance(grouped_sample_hamming_by_length, dict):
        raise ValueError(f"{grouped_sample_hamming_by_length} must be a dictionary.")

    min_distances = {}

    for length, comparisons in grouped_sample_hamming_by_length.items():
        distances = [dist for _, _, dist in comparisons]
        if distances:
            min_distances[length] = min(distances)

    return min_distances


def assert_min_hamming_above_threshold(min_distances_by_group: dict, max_hamming: int) -> None:
    """
    Validates that the minimum Hamming distance in each group is >= max_hamming.

    Args:
        min_distances_by_group (dict): Output from `find_min_hamming_distances`,
            mapping group length to its minimum Hamming distance.
        max_hamming (int): The minimum allowed Hamming distance.

    Raises:
        ValueError: If any group has a minimum Hamming distance less than max_hamming.
    """
    if not isinstance(min_distances_by_group, dict):
        raise ValueError(f"{min_distances_by_group} must be a dictionary.")
    if not isinstance(max_hamming, int):
        raise ValueError(f"{max_hamming} must be an integer.")

    for length, min_distance in min_distances_by_group.items():
        if min_distance < max_hamming:
            raise ValueError(f"Minimum Hamming distance violation in group {length}: " f"{min_distance} < {max_hamming}")


def custom_priority_by_length_sort_key(key):
    """
    Sorts tuple or integer keys based on custom priority: keys without zeros are prioritized
    and sorted by descending total length, while keys containing a zero are placed last.

    This method is intended for sorting barcode length indicators. If a single integer is
    passed instead of a tuple, it is treated as (value, 0).

    Args:
        key (tuple or int): A tuple of two integers representing barcode segment lengths,
            or a single integer to be treated as (value, 0).

    Returns:
        tuple: A tuple used for sorting. The first element indicates whether the key
            contains a zero (False = higher priority), and the second is the negative
            sum of the tuple values (for descending order).

    Raises:
        TypeError: If the input is not an int or a tuple of one or two integers.
    """
    if key is None:
        raise TypeError("Key cannot be None.")

    if isinstance(key, int):
        key = (key, 0)
    elif isinstance(key, tuple):
        if len(key) == 1:
            key = (key[0], 0)
        elif len(key) != 2:
            raise TypeError("Tuple key must have one or two integers.")
    else:
        raise TypeError("Input must be an int or a tuple of one or two integers.")

    if not all(isinstance(i, int) for i in key):
        raise TypeError("All elements in key must be integers.")

    has_zero = 0 in key
    total_len = sum(key)
    return (has_zero, -total_len)


def extract_index_from_header_illumina(name: str) -> str:
    """
    Extract the index sequence from a FASTQ read header.

    The function assumes that the index is located at the end of the header string,
    separated by spaces and colons. It extracts the portion of the header after the
    third colon in the last space-separated segment.

    Args:
        name (str): The read header from the FASTQ file. This is typically a string
                    containing metadata about the read, including the index.

    Returns:
        str: The extracted index sequence from the header.

    Raises:
        ValueError: If the header is empty, None, or malformed (e.g., contains too
                    many colons in the index portion).
    """
    if not name:
        raise ValueError("Read header is empty.")
    if name is None:
        raise ValueError("Read header is None.")

    # Split the name by spaces and take the last part
    # this assumes the index is always at the end of the header
    split_str = name.rsplit(" ", 1)[-1]
    # Split that last part by colons and join everything after the third colon
    parts = split_str.split(":")
    index_part = ":".join(parts[3:])

    # Check for malformed index
    if index_part.count(":") > 1:
        raise ValueError("Too many colons in index portion — possibly malformed index.")

    return index_part


def trim_merge_string(input_str: str, length: int) -> str:
    """
    Trims and optionally splits and merges a string based on a specified length.

    If the input string contains no space, it is trimmed from the end to the specified length.
    If the input string contains a space, it is split into two parts at the first space. Each part
    is trimmed to `length // 2` characters, and the two trimmed parts are concatenated.

    Args:
        length (int): The total number of characters to retain.
        input_str (str): The input string to process.

    Returns:
        str: A trimmed string (if no space is present), or a merged string composed of two trimmed
            parts (if a space is present).

    Raises:
        ValueError: If the provided length is negative.
    """
    if input_str is None:
        raise ValueError("Input string is None.")
    if length is None:
        raise ValueError("Length is None.")

    if not isinstance(input_str, str):
        raise ValueError(f"{input_str} must be a string.")
    if not isinstance(length, int):
        raise ValueError(f"{length} must be an integer.")

    if length < 0:
        raise ValueError("Length must be non-negative.")

    # check if the string contains a space
    if " " in input_str:
        part1, part2 = input_str.split(" ", 1)
        half_len = length // 2  # Always returns a float
        return part1[:half_len] + part2[:half_len]
    else:
        return input_str[:length]


def build_bk_tree_index(sorted_grouped_samples_by_length: dict, grouped_samples_by_length: dict) -> dict:
    """
    Builds a nested dictionary of BK-trees for index matching based on index length and sample identifiers.

    This method constructs BK-tree data structures for each sample's index information, grouped by index length.
    It iterates over the provided index lengths and corresponding sample data, creating BK-trees for the first and
    optional second index sequences per sample. The resulting structure is organized by index length, then sample ID,
    with each entry containing a list of one or two BK-trees.

    Args:
        sorted_grouped_samples_by_length (dict): A dictionary where keys are index lengths (int) and values are
            ordered lists or other iterable forms used to determine processing order.
        grouped_samples_by_length (dict): A nested dictionary of the form
            {length: {sample_id: [index1_seq, index2_seq (optional)]}}, where each index sequence is a list of strings.

    Returns:
        dict: A nested dictionary of the form
            {length: {sample_id: [BKTree(index1), BKTree(index2 or None)]}}, used for approximate matching
            of barcodes using Hamming distance.

    Raises:
        TypeError: If input structures are not correctly formatted or contain unexpected data types.
    """
    grouped_bk_trees = {}
    for length in sorted_grouped_samples_by_length:
        if length not in grouped_bk_trees:
            grouped_bk_trees[length] = {}
        for sample in grouped_samples_by_length[length]:
            index1_tree = BKTree(hamming_distance, [grouped_samples_by_length[length][sample][0]])
            index2_tree = (
                BKTree(hamming_distance, [grouped_samples_by_length[length][sample][1]])
                if len(grouped_samples_by_length[length][sample]) > 1
                else None
            )
            grouped_bk_trees[length][sample] = [index1_tree, index2_tree]

    return grouped_bk_trees


# def index_to_match_key(read_header: str, barcode_bktree_map: dict, max_hamming_distance_1: int, max_hamming_distance_2: int = None) -> tuple:
#     """
#     Matches a read index to a sample using BK-tree search.

#     Args:
#         read_header (str): Read name/header containing the barcode.
#         sorted_group_lengths (list): Sorted list of index length tuples.
#         grouped_bk_trees (dict): Maps (i5_len, i7_len) to BK-trees and barcode mapping.
#         max_hamming_distance (int): Max allowable per-index Hamming distance.

#     Returns:
#         tuple: (matched sample name or 'undetermined', trimmed index used)
# """
# raw_index = extract_index_from_header_illumina(read_header)
# index_split = re.split(r"[^A-Za-z]", raw_index)
# read_index1, read_index2 = (index_split[0], index_split[1]) if len(index_split) > 1 else (index_split[0], None)

# best_index1_match = None
# best_index2_match = None

# index1_exact_match_found = False
# index2_exact_match_found = False

# sample_index_match = {}
# # best_total = max_hamming_distance_1 + 1  # Initialize to a value greater than max_hamming_distance
# # best_sample = "undetermined"

# for length_key in barcode_bktree_map:
#     read_index1_trimmed = trim_merge_string(read_index1, length_key[0])
#     read_index2_trimmed = trim_merge_string(read_index2, length_key[1]) if length_key[1] > 0 else None

#     for sample, (index1_tree, index2_tree) in barcode_bktree_map[length_key].items():
#         # Only search if we haven't found a perfect match already

#         # best_index1 = 0
#         # best_index2 = 0

#         if not index1_exact_match_found:
#             index1_matches = index1_tree.find(read_index1_trimmed, max_hamming_distance_1)
#     for dist, val in index1_matches:
#         #     print(dist, val)
#         # if dist == 0:
#         #     best_index1_match = [val, sample]
#         #     index1_exact_match_found = True
#         # break  # stop checking index1 matches
#         # if best_index1_match is None or dist < best_index1_match[0]:
#         if best_index1_match is None:
#             best_index1_match = [dist, sample]
#         if best_index1_match is not None and dist < best_index1_match[0]:
#             best_index1_match = [dist, sample]

#     # best_index1 = min(index1_matches, key=lambda x: x[0])

# # print("best_index1_match: ", best_index1_match)

# if not index2_exact_match_found and length_key[1] > 0:
#     if max_hamming_distance_2 is None:
#         index2_matches = index2_tree.find(read_index2_trimmed, max_hamming_distance_1)
#     else:
#         index2_matches = index2_tree.find(read_index2_trimmed, max_hamming_distance_2)
#     for dist2, val2 in index2_matches:
#         print(dist2)
#         if dist2 == 0:
#             best_index2_match = [dist2, sample]
#             index2_exact_match_found = True
#             break  # stop checking index2 matches
#         if best_index2_match is None or dist2 < best_index2_match[0]:
#             best_index2_match = [dist2, sample]
# print("best_index2_match: ", best_index2_match)

# sample_index_match[sample] = (best_index1_match, best_index2_match)

# best_index2 = min(index2_matches, key=lambda x: x[0]) if index2_matches else (max_hamming_distance_2 + 1, None)

#     total_distance = best_index1[0] + best_index2[0]
#     if total_distance < best_total:
#         best_total = total_distance
#         best_sample = sample

#     print(total_distance, best_index1, best_index2)

# return best_sample

#         # Break the sample loop if exact matches for both indexes are found
#         if index1_exact_match_found and index2_exact_match_found:
#             break

#     # Break the barcode_bktree_map loop if exact matches for both indexes are found
#     if index1_exact_match_found and index2_exact_match_found:
#         break

# # If no matches were found, return "undetermined"
# if best_index1_match is None and best_index2_match is None:
#     return "undetermined"

# print("here: " + str(best_index1_match) + " " + str(best_index2_match))
# # If both indexes have been matched to the same sample, return the sample name
# if best_index1_match[1] == best_index2_match[1]:
#     print("Returning matched sample: " + best_index1_match[1])
#     return best_index1_match[1]  # Return sample name
# print(f"Trying to match i5: {index1_matches}, i7: {index2_matches}  against length_key: {length_key}")

# return best_index1_match, best_index2_match


# def index_to_match_key(read_header: str, barcode_bktree_map: dict, max_hamming_distance_1: int, max_hamming_distance_2: int = None) -> tuple:
#     """
#     Matches a read index to a sample using BK-tree search.

#     Args:
#         read_header (str): Read name/header containing the barcode.
#         sorted_group_lengths (list): Sorted list of index length tuples.
#         grouped_bk_trees (dict): Maps (i5_len, i7_len) to BK-trees and barcode mapping.
#         max_hamming_distance (int): Max allowable per-index Hamming distance.

#     Returns:
#         tuple: (matched sample name or 'undetermined', trimmed index used)
#     """
#     raw_index = extract_index_from_header_illumina(read_header)
#     index_split = re.split(r"[^A-Za-z]", raw_index)
#     read_index1, read_index2 = (index_split[0], index_split[1]) if len(index_split) > 1 else (index_split[0], None)

#     # # If we have a separate hamming distance value for each index, the total combined distance is the sum of both, else we use the only hamming distance value provided as the default value for both indexes
#     # max_hamming_distance_combined = max_hamming_distance_1 + (max_hamming_distance_2 if max_hamming_distance_2 is not None else max_hamming_distance_1) + 1
#     # # any possible matches must have a hamming distance equal to or less than (max_hamming_distance_1+max_hamming_distance_2), so we initialize the best matches to a value greater than that
#     # best_index1_match = [max_hamming_distance_combined, "undetermined"]
#     # best_index2_match = [max_hamming_distance_combined, "undetermined"]

#     best_index1_match = [max_hamming_distance_1 + 1, "undetermined"]
#     best_index2_match = [max_hamming_distance_2 + 1, "undetermined"] if max_hamming_distance_2 else [max_hamming_distance_1 + 1, "undetermined"]

#     index1_exact_match_found = False
#     index2_exact_match_found = False

#     sample_index_match = {}

#     for length_key in barcode_bktree_map:
#         read_index1_trimmed = trim_merge_string(read_index1, length_key[0])
#         read_index2_trimmed = trim_merge_string(read_index2, length_key[1]) if length_key[1] > 0 else None

#         for sample, (index1_tree, index2_tree) in barcode_bktree_map[length_key].items():
#             # Only search if we haven't found a perfect match already
#             if not index1_exact_match_found:
#                 index1_matches = index1_tree.find(read_index1_trimmed, max_hamming_distance_1)
#                 for dist, val in index1_matches:
#                     if dist < best_index1_match[0]:
#                         best_index1_match = [dist, sample]

#             if not index2_exact_match_found and length_key[1] > 0:
#                 if max_hamming_distance_2 is None:
#                     index2_matches = index2_tree.find(read_index2_trimmed, max_hamming_distance_1)
#                 else:
#                     index2_matches = index2_tree.find(read_index2_trimmed, max_hamming_distance_2)
#                 for dist2, val2 in index2_matches:
#                     # print(dist2)
#                     if dist2 == 0:
#                         best_index2_match = [dist2, sample]
#                         index2_exact_match_found = True
#                         break  # stop checking index2 matches
#                     if best_index2_match[1] == "undetermined" or dist2 < best_index2_match[0]:
#                         best_index2_match = [dist2, sample]

#             # Break the sample loop if exact matches for both indexes are found
#             if index1_exact_match_found and index2_exact_match_found:
#                 break

#             if sample == "sample_4" and best_index1_match[0] < 2 and best_index2_match[0] < 2:
#         # #     print(f"Trying to match i5: {index1_matches}, i7: {index2_matches}  against length_key: {length_key}")
#                 print("best_index1_match: ", best_index1_match, best_index2_match, sample)

#         # Break the barcode_bktree_map loop if exact matches for both indexes are found
#         if index1_exact_match_found and index2_exact_match_found:
#             break

#     if read_index2 is None:
#         return best_index1_match[1]
#     else:
#         # If no matches were found, return "undetermined"
#         if best_index1_match[1] == "undetermined" and best_index2_match[1] == "undetermined" or best_index1_match[1] != best_index2_match[1]:
#             return "undetermined"

#         # If both indexes have been matched to the same sample, return the sample name
#         if best_index1_match[1] == best_index2_match[1]:
#             if best_index1_match[1] == "sample_4":
#                 print("Returning matched sample: " + best_index1_match[1])
#             return best_index1_match[1]  # Return sample name


def index_to_match_key(read_header: str, barcode_bktree_map: dict, max_hamming_distance_1: int, max_hamming_distance_2: int = None) -> str:
    """
    Matches a read index to a sample using BK-tree search based on combined Hamming distances.

    This function extracts index barcodes from the read header, searches them against a BK-tree
    structure using configurable Hamming distance thresholds, and identifies the most likely
    sample match based on minimal combined distance.

    Args:
        read_header (str): Read name/header containing the barcode.
        barcode_bktree_map (dict): Maps (i5_len, i7_len) to BK-trees and sample barcode mappings.
        max_hamming_distance_1 (int): Maximum Hamming distance allowed for the first index (i5).
        max_hamming_distance_2 (int, optional): Maximum Hamming distance for the second index (i7).
            Defaults to the value of `max_hamming_distance_1` if not provided.

    Returns:
        str: Matched sample name, or 'undetermined' if no valid match found.
    """
    raw_index = extract_index_from_header_illumina(read_header)
    index_split = re.split(r"[^A-Za-z]", raw_index)
    read_index1, read_index2 = (index_split[0], index_split[1]) if len(index_split) > 1 else (index_split[0], None)

    if max_hamming_distance_2 is None:
        max_hamming_distance_2 = max_hamming_distance_1

    sample_match_scores = {}

    for length_key in barcode_bktree_map:
        read_index1_trimmed = trim_merge_string(read_index1, length_key[0])
        read_index2_trimmed = trim_merge_string(read_index2, length_key[1]) if read_index2 and length_key[1] > 0 else None

        for sample, (index1_tree, index2_tree) in barcode_bktree_map[length_key].items():
            index1_dist = None
            index2_dist = None

            # Search index1
            index1_matches = index1_tree.find(read_index1_trimmed, max_hamming_distance_1)
            for dist, _ in index1_matches:
                index1_dist = dist
                break  # Take the closest match

            # Search index2 if applicable
            if read_index2_trimmed and length_key[1] > 0:
                index2_matches = index2_tree.find(read_index2_trimmed, max_hamming_distance_2)
                for dist2, _ in index2_matches:
                    index2_dist = dist2
                    break  # Take the closest match

            # If at least one index matched
            if index1_dist is not None or index2_dist is not None:
                index1_dist = index1_dist if index1_dist is not None else max_hamming_distance_1 + 1
                index2_dist = index2_dist if index2_dist is not None else max_hamming_distance_2 + 1

                if index1_dist <= max_hamming_distance_1 and index2_dist <= max_hamming_distance_2:
                    sample_match_scores[sample] = index1_dist + index2_dist

    if not sample_match_scores:
        return "undetermined"

    # Return sample with the lowest total Hamming distance
    best_sample = min(sample_match_scores.items(), key=lambda x: x[1])[0]

    return best_sample


# def find_closest_match(barcode_dict: dict, seq: str, max_hamming: int) -> str:
#     """
#     Finds the sample and barcode with the smallest Hamming distance to the given sequence.

#     Args:
#         barcode_dict (dict): Dictionary mapping sample names to barcode strings.
#         seq (str): The sequence to compare against the barcodes.
#         max_hamming (int): Maximum allowed Hamming distance.

# Returns:
#     str or None: Returns sample_name with the smallest
#         Hamming distance, or None if no barcode is within max_hamming.
# """
# if not isinstance(seq, str):
#     raise ValueError(f"{seq} must be a string.")

# if not isinstance(barcode_dict, dict):
#     raise ValueError(f"{barcode_dict} must be a dictionary.")

# if not isinstance(max_hamming, int):
#     raise ValueError(f"{max_hamming} must be an integer.")

# if max_hamming < 0:
#     raise ValueError(f"{max_hamming} must be a non-negative integer.")

# best_match = "undetermined"
# min_distance = float("inf")

# for sample, barcode in barcode_dict.items():
#     dist = hamming_distance(seq, barcode)
#     # If the strings are a perfect match, skip further comparisons
#     if dist == 0:
#         return sample

#     if dist < min_distance and dist <= max_hamming:
#         min_distance = dist
#         best_match = sample

# return best_match


# def build_bk_trees(barcode_sample_dict):
#     grouped = defaultdict(list)
#     for sample, barcode in barcode_sample_dict.items():
#         for bc in barcode.split(','):
#             grouped[len(bc)].append((bc, sample))  # (barcode, sample)
#         print(grouped)

# trees = {}
# for length, entries in grouped.items():
#     tree = BKTree(hamming_distance)
#     for bc, _ in entries:
#         tree.add(bc)
#     trees[length] = (tree, {bc: sample for bc, sample in entries})
# return trees


# def find_matching_sample(barcode, trees, max_hamming_distance):
#     length = len(barcode)
#     if length not in trees:
#         return None
#     tree, barcode_to_sample = trees[length]
#     matches = tree.find(barcode, max_hamming_distance)
#     if matches:
#         # Prefer lowest distance match
#         matches.sort()
#         matched_barcode = matches[0][1]
#         return barcode_to_sample[matched_barcode]
#     return None


def demultiplex_fastq_by_barcode(
    samples_barcode_from_dict: dict, fastq_file_r1: str, max_hamming_distance: int = 0, output_dir: str = ".", fastq_file_r2: str = None
) -> dict:
    """
    Demultiplexes a FASTQ file by assigning reads to samples based on barcode sequences.

    This function takes a FASTQ file and a dictionary of sample barcodes, groups the samples
    by the length of their barcodes, and then assigns reads to samples based on the closest
    barcode match within a specified Hamming distance threshold. Reads that do not match any
    barcode within the threshold are written to an "undetermined" file.

    The method enforces barcode dissimilarity by checking that the minimum Hamming distance
    between barcodes in each group exceeds the specified threshold.

    Args:
        fastq_file (str): Path to the FASTQ file to be demultiplexed.
        samples_barcode_from_dict (dict): A dictionary mapping sample names to their barcodes,
            which may be simple strings or nested dicts with keys like "index" or "index2".
        max_hamming_distance (int, optional): Maximum allowable Hamming distance between
            a read’s index and a sample barcode to be considered a match. Defaults to 0.
        output_dir (str, optional): Directory where demultiplexed FASTQ files will be written.
            Defaults to the current directory.

    Returns:
        dict: A dictionary mapping sample names to the number of reads assigned to each.

    Raises:
        ValueError: If any group of barcodes contains pairs with a minimum Hamming distance
            less than or equal to the provided `max_hamming_distance`.
        TypeError: If the barcode structure is not a valid string or expected dict format.
        FileNotFoundError: If the input FASTQ file does not exist or cannot be opened.
        IOError: If any of the output files cannot be created or written to.
    """
    ## Group samples by index length
    grouped_samples_by_length = group_samples_by_index_length(samples_barcode_from_dict)

    # Compare the all the barcodes in each group against each other to find the Hamming distance for each pair compared
    grouped_sample_by_length_hamming_value = {}
    grouped_sample_by_length_hamming_value = {
        length: crosscheck_barcode_proximity(samples)
        for length, samples in grouped_samples_by_length.items()
        if len(samples) > 1  # Skip groups with only one entry
    }
    # print(grouped_sample_by_length_hamming_value)
    # find the minimum hamming distance for each index-length group
    min_hamming_distances_by_length = find_min_hamming_distances(grouped_sample_by_length_hamming_value)
    # check if the minimum hamming distance is above the threshold max hamming distance
    # returns ValueError if any group has a minimum Hamming distance less than max_hamming
    assert_min_hamming_above_threshold(min_hamming_distances_by_length, max_hamming_distance)
    if max_hamming_distance_2 is not None:
        assert_min_hamming_above_threshold(min_hamming_distances_by_length, max_hamming_distance_2)
    # print(assert_min_hamming_above_threshold(min_hamming_distances_by_length, max_hamming_distance))

    ## Sort the grouped samples by length, prioritizing those without zeros
    sorted_grouped_samples_by_length = sorted(grouped_samples_by_length.keys(), key=custom_priority_by_length_sort_key)
    # print(sorted_grouped_samples_by_length)

    ## Create a fastq file for each sample + an "undetermined" file for unassigned reads
    if isinstance(fastq_file_r1, str):
        # Extract read and lane information from the R1 fastq file name
        read_lane_r1_info = re.search(r"S\d+_([^\.]+)", fastq_file_r1)
        if not read_lane_r1_info:
            read_lane_r1_info = fastq_file_r1.split("_", 1)[1].split(".", 1)[0]
    file_handles_r1 = {
        sample: gzip.open(os.path.join(output_dir, f"{sample}_{read_lane_r1_info}.fastq.gz"), "ab") for sample in samples_barcode_from_dict
    }
    file_handles_r1["undetermined"] = gzip.open(os.path.join(output_dir, f"undetermined_{read_lane_r1_info}.fastq.gz"), "ab")

    if fastq_file_r2:
        if isinstance(fastq_file_r2, str):
            # Extract read and lane information from the R2 fastq file name
            read_lane_r2_info = re.search(r"S\d+_([^\.]+)", fastq_file_r2)
            if not read_lane_r2_info:
                read_lane_r2_info = fastq_file_r2.split("_", 1)[1].split(".", 1)[0]
        file_handles_r2 = {
            sample: gzip.open(os.path.join(output_dir, f"{sample}_{read_lane_r2_info}.fastq.gz"), "ab") for sample in samples_barcode_from_dict
        }
        file_handles_r2["undetermined"] = gzip.open(os.path.join(output_dir, f"undetermined_{read_lane_r2_info}.fastq.gz"), "ab")

    # Keep track of each read assigned to each sample
    sample_assigned_read = defaultdict(list)  # this one keeps track of sample+read index information
    sample_count = {}  # this one keeps track of how many reads are assigned to each sample

    ## Read the FASTQ file and extract indexes, save all indexes in a list
    fastq_1 = FastqFile(fastq_file_r1)
    if fastq_file_r2:
        fastq_2 = FastqFile(fastq_file_r2)
        fastq_2_iter = fastq_2.open_read_iterator(as_string=True)

    # Build BK-trees based on a pre-determined Hamming sequence
    # barcode_map = {}
    grouped_bk_trees = build_bk_tree_index(sorted_grouped_samples_by_length, grouped_samples_by_length)
    # for length in sorted_grouped_samples_by_length:
    #     if length not in barcode_map:
    #         barcode_map[length] = {}
    #     # print(f"Processing length group: {grouped_samples_by_length[length]}")
    #     for sample in grouped_samples_by_length[length]:
    #         # print(f"Building BK-tree for length {length} with samples: {grouped_samples_by_length[length]}")
    #         # print(grouped_samples_by_length[length][sample])
    #         # index1_tree, index2_tree, barcode_map = BKTree(hamming_distance, grouped_samples_by_length[length][sample])
    #         index1_tree = BKTree(hamming_distance, grouped_samples_by_length[length][sample][0])
    #         index2_tree = BKTree(hamming_distance, grouped_samples_by_length[length][sample][1]) if len(grouped_samples_by_length[length][sample]) > 1 else None
    #         # print(index2_tree)
    #         barcode_map[length][sample] = [index1_tree, index2_tree]

    # print(grouped_samples_by_length)

    for name, seq, qual in fastq_1.open_read_iterator(as_string=True):
        # print(f"name: {name}")
        if max_hamming_distance_2 is not None:
            # Use the index_to_match_key function to find the best match for the read
            match = index_to_match_key(name, grouped_bk_trees, max_hamming_distance, max_hamming_distance_2)
        else:
            match = index_to_match_key(name, grouped_bk_trees, max_hamming_distance)
        # print(match)

        # raw_index = extract_index_from_header_illumina(name)
        #     match = find_matching_sample(raw_index, trees, max_hamming_distance) or 'undetermined'

        # Assign the read to the matched sample
        sample_assigned_read[match].append(name)

        # Write the sequence to the appropriate file
        FastqFile.write_read(file_handles_r1[match], name, seq, qual)

        if fastq_file_r2:
            # Read the corresponding R2 read
            name_r2, seq_r2, qual_r2 = next(fastq_2_iter)
            # Write the R2 read to the appropriate file
            FastqFile.write_read(file_handles_r2[match], name_r2, seq_r2, qual_r2)

    for sample, reads in sample_assigned_read.items():
        sample_count[sample] = len(reads)

    for f in file_handles_r1.values():
        f.close()

    if fastq_file_r2:
        for f in file_handles_r2.values():
            f.close()

    return sample_count
