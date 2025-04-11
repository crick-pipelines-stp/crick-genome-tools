import itertools
import re

import numpy as np

from crick_genome_tools.io.fastq_file import FastqFile


# from carmack.barcode.barcode_extractor import BarcodeExtractor


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


def group_samples_by_index_length(sample_index_dict: dict) -> dict:
    """
    Groups samples by the length of their barcode indices, counting only alphabetic characters.

    This function takes a dictionary of samples and their associated barcodes or nested metadata.
    If a sample's value is a string, it's treated directly as a barcode. If it's a dictionary,
    the function searches for a "barcode" or "index" key to extract the barcode. The barcodes are
    then grouped into a new dictionary based on the number of alphabetic characters they contain.

    Args:
        sample_index_dict (dict): A dictionary where each key is a sample name (str) and each value
            is either a barcode string or a dictionary containing at least a "barcode" or "index" key.

    Returns:
        dict: A dictionary where each key is an integer representing the count of alphabetic characters
            in the barcode, and the corresponding value is a dictionary of sample-barcode pairs that
            share that alphabetic character count.

    Raises:
        TypeError: If the input is not a dictionary, or if any barcode is not a string.
        KeyError: If a sample has a dict value but lacks both "barcode" and "index" keys.
    """
    if not isinstance(sample_index_dict, dict):
        raise TypeError("Input must be a dictionary.")

    grouped = {}
    for sample, value in sample_index_dict.items():
        if isinstance(value, str):
            # if sample has a string value, treat it as a barcode
            barcode = value

        elif isinstance(value, dict):
            # if sample has a dict value, look for "barcode" or "index" keys
            if "barcode" in value:
                barcode = value["barcode"]
            elif "index" in value:
                barcode = value["index"]
            else:
                barcode = ""
        else:
            raise TypeError(f"Sample '{sample}' has an unsupported value type: {type(value).__name__}")

        if not isinstance(barcode, str):
            raise TypeError(f"Barcode for sample '{sample}' must be a string.")

        # Calculate the length of the barcode, while ignoring non-alphabetic characters
        barcode_len = sum(c.isalpha() for c in barcode)
        grouped.setdefault(barcode_len, {})[sample] = barcode

    return grouped


##################################################################################
DNA_ALPHABET = "AGCT"
ALPHABET_MINUS = {char: {c for c in DNA_ALPHABET if c != char} for char in DNA_ALPHABET}  # This is a set of alternative bases given a base
ALPHABET_MINUS["N"] = set(DNA_ALPHABET)


def gen_nearby_seqs(seq, barcode_set, maxdist):
    """Generate all sequences with at most maxdist changes from seq that are in a provided seq, along with the
    quality values of the bases at the changed positions. Automatically will target N's in a sequence as letters
    which must be changed. If there are more N's than allowed changes - we return nothing
    """

    new_seq = set()

    # Find all index positions which are not N in seq as a list
    non_n_indices = [i for i in range(len(seq)) if seq[i] != "N"]

    # Find all positions which are N in seq as a tuple
    n_indices = tuple([i for i in range(len(seq)) if seq[i] == "N"])

    # The number of unknown N's dicates the minmimum hamming distance that combinations must be from the original sequence
    mindist = len(n_indices)

    # If this is too far away then we just return nothing
    if mindist > maxdist:
        return [], 0

    # If the input sequence is in the barcode set, include the seq and qs in the output
    if seq in barcode_set:
        yield seq

    # Combinations are generated in batches by changing n number of indices in the sequence, then n+1 and so on
    # The min number of positions to change is dictated by the number of N's in the sequence
    # The max number of positions to change is dictated by the max hamming distance
    for dist in range(mindist, maxdist + 1):

        # Generate possible combinations of non-required indices to change for this hamming distance level
        # This list will be empty if the number of N's is equal to the hamming distance
        for modified_indices in itertools.combinations(non_n_indices, dist - mindist):

            # Combine the indices we have to change because of N's and the other potential cominations into a final
            # List of indices to change
            indices = set(modified_indices + n_indices)

            # Convert the set to a list of indices for subsetting the qs scores (ignore the empty list at the beggining)
            indices_list = np.array(list(indices))
            if len(indices_list) == 0:
                continue

            # Generate possible base substitutions from the indice positions using the minus alphabet
            for substitutions in itertools.product(*[ALPHABET_MINUS[base] if i in indices else base for i, base in enumerate(seq)]):
                new_seq = "".join(substitutions)

                # If the new sequence is in the whitelist, sum the QS scores for the changed sequences and return
                yield new_seq


##################################################################


def generate_nearby_barcodes_by_length(grouped_barcodes: dict, max_hamming_distance: str) -> dict:
    """
    Generates nearby barcodes within a specified Hamming distance for each sample in grouped barcode data.

    This function processes a dictionary where keys represent barcode lengths and values are dictionaries
    mapping sample names to barcode sequences. It uses the `gen_nearby_seqs` function to identify all
    barcode sequences that differ from the original by at most a given Hamming distance and are present
    within the same barcode length group.

    Dual-index barcodes (e.g., "ACGT+AGGT") are detected automatically and returned as
    concatenated strings (e.g., "ACGTAGGT") with all valid nearby combinations.

    Args:
        grouped_barcodes (dict): Dictionary with barcode length as keys and inner dicts of
                                 sample names to barcodes (single or dual-index).
        max_hamming_distance (int): Maximum allowed Hamming distance.

    Returns:
        dict: {length: {sample: set of nearby full-barcode strings}}

    Raises:
        TypeError, ValueError: If inputs are malformed.
    """
    if not isinstance(grouped_barcodes, dict):
        raise TypeError("Input must be a dictionary.")

    if not isinstance(max_hamming_distance, int) or max_hamming_distance < 0:
        raise ValueError("max_hamming_distance must be a non-negative integer.")

    result = {}

    for length, samples in grouped_barcodes.items():
        if not isinstance(samples, dict):
            raise TypeError(f"{grouped_barcodes} dict is incorrectly formatted.")

        # Build a set of all barcode components
        barcode_set = set()
        for barcode in samples.values():
            parts = [p for p in re.split(r"[^A-Za-z]+", barcode) if p]
            barcode_set.update(parts)

        near_matches = {}

        for sample, barcode in samples.items():
            parts = [p for p in re.split(r"[^A-Za-z]+", barcode) if p]

            # Generate nearby variants for each part
            matches_sets = []
            for part in parts:
                matches = set()
                for near_seq in gen_nearby_seqs(part, barcode_set, max_hamming_distance):
                    matches.add(near_seq)
                matches_sets.append(matches)

            # Combine parts back into full barcodes (as one string)
            if len(matches_sets) > 1:
                # Dual-index: combine cross-product of variants
                combined = {"".join(combo) for combo in itertools.product(*matches_sets)}
            else:
                # Single-index
                combined = matches_sets[0]

            near_matches[sample] = combined

        result[length] = near_matches

    return result


# def read_fastq(self, fastq_file):
#     """
#     Read a FASTQ file and return a list of sequences.

#     Args:
#         file_path (str): Path to the FASTQ file.

#     Returns:
#         list: List of sequences.
#     """
#     fastq = FastqFile(fastq_file)

#     for name, seq, qual in fastq.open_read_iterator(as_string=True):
#         self.extract_index_from_header(name)


def gen_index_and_hamming_dict(fastq_file, group1, group2):
    """
    Generate index and Hamming distance for a FASTQ file.

    Args:
        fastq_file (str): Path to the FASTQ file.
        group1 (dict): Sample names mapped to lists of barcodes.
        group2 (list): List of dicts, each with 'barcode' and 'data' keys.

    Returns:
        None
    """
    fastq = FastqFile(fastq_file)

    for name, seq, qual in fastq.open_read_iterator(as_string=True):
        index = extract_index_from_header_illumina(name)
        print(f"Index: {index}, Sequence: {seq}, Quality: {qual}")


#############################################################################

import os
from collections import defaultdict


def build_barcode_maps(group1):
    """
    Builds two maps:
        - exact_lookup: barcode → (sample, is_dual, length)
        - partial_info: list of (sample, partials, is_dual, length)

    Args:
        group1 (dict): Sample names mapped to lists of barcodes.

    Returns:
        tuple: (exact_lookup dict, partial_info list)
    """
    exact_lookup = {}
    partial_info = []

    for sample, barcodes in group1.items():
        for bc in barcodes:
            is_dual = "-" in bc
            parts = bc.split("-")
            length = sum(len(p) for p in parts)
            exact_lookup[bc] = (sample, is_dual, length)
            partial_info.append((sample, parts, is_dual, length))

    return exact_lookup, partial_info


def score_match(is_exact, is_dual, length):
    """
    Scores a match based on exactness and complexity.

    Args:
        is_exact (bool): True if exact match.
        is_dual (bool): True if dual-indexed barcode.
        length (int): Total length of the barcode.

    Returns:
        int: Match score.
    """
    score = 0
    if is_exact:
        score += 1000
    if is_dual:
        score += 100
    score += length
    return score


def compare_and_write(group1, group2, output_dir):
    """
    Compares reads to sample barcodes using exact-first, complexity-aware matching.
    Writes matched reads to sample files. Unmatched reads go to "undetermined.txt".

    Args:
        group1 (dict): Sample names mapped to lists of barcodes.
        group2 (list): List of dicts, each with 'barcode' and 'data' keys.
        output_dir (str): Directory to write output files.

    Returns:
        None

    Raises:
        FileNotFoundError: If output directory cannot be created.
    """
    os.makedirs(output_dir, exist_ok=True)

    exact_lookup, partial_info = build_barcode_maps(group1)

    file_handles = defaultdict(lambda: open(os.path.join(output_dir, "undetermined.txt"), "a"))
    for sample in group1:
        file_handles[sample] = open(os.path.join(output_dir, f"{sample}.txt"), "a")

    for read in group2:
        read_bc = read["barcode"]
        best_sample = "undetermined"
        best_score = -1

        # ✅ Step 1: Exact match
        if read_bc in exact_lookup:
            sample, is_dual, length = exact_lookup[read_bc]
            best_sample = sample
            best_score = score_match(True, is_dual, length)

        else:
            # 🔁 Step 2: Partial matches
            for sample, parts, is_dual, length in partial_info:
                if read_bc in parts:  # partial match to index1 or index2
                    score = score_match(False, is_dual, length)
                    if score > best_score:
                        best_score = score
                        best_sample = sample

        file_handles[best_sample].write(read["data"] + "\n")

    # ✅ Close all files
    for f in file_handles.values():
        f.close()
