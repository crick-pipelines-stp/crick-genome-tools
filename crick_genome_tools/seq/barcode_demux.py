import itertools
import re

import numpy as np

from crick_genome_tools.io.fastq_file import FastqFile


DNA_ALPHABET = "AGCT"
ALPHABET_MINUS = {char: {c for c in DNA_ALPHABET if c != char} for char in DNA_ALPHABET}  # This is a set of alternative bases given a base
ALPHABET_MINUS["N"] = set(DNA_ALPHABET)


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


def gen_nearby_seqs(seq: str, barcode_set, maxdist: int = 0) -> str:
    """
    Generate all sequences with at most maxdist changes from seq that are in a provided seq, along with the
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
        return []

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
            for substitutions in itertools.product(*[ALPHABET_MINUS[base.upper()] if i in indices else base for i, base in enumerate(seq)]):
                new_seq = "".join(sub.upper() if original_base.isupper() else sub.lower() for sub, original_base in zip(substitutions, seq))
                # If the new sequence is in the whitelist, sum the QS scores for the changed sequences and return
                yield new_seq

def generate_nearby_barcodes_by_length(grouped_barcodes: dict, max_hamming_distance: str) -> dict:
    """
    Generates nearby barcodes within a specified Hamming distance for each sample in grouped barcode data.

    This function processes a dictionary where keys represent barcode lengths and values are dictionaries
    mapping sample names to barcode sequences. It uses the `gen_nearby_seqs` function to identify all
    barcode sequences that differ from the original by at most a given Hamming distance and are present
    within the same barcode length group.

    Treats dual indexes (e.g. 'ACGT+AGGT') as single merged barcodes (e.g. 'ACGTAGGT').

    Args:
        grouped_barcodes (dict): Nested dict {length: {sample: barcode or dual-barcode string}}.
        max_hamming_distance (int): Maximum number of allowed mutations.

    Returns:
        dict: Nested dict {length: {sample: set of nearby barcode strings}}.

    Raises:
        TypeError: If inputs are not properly formatted.
        ValueError: If max_hamming_distance is invalid.
    """

    if not isinstance(grouped_barcodes, dict):
        raise TypeError("Input must be a dictionary.")

    if not isinstance(max_hamming_distance, int) or max_hamming_distance < 0:
        raise ValueError("max_hamming_distance must be a non-negative integer.")

    result = {}

    for length, samples in grouped_barcodes.items():
        if not isinstance(samples, dict):
            raise TypeError(f"{grouped_barcodes} dict is incorrectly formatted.")

        # Remove separators from dual index barcodes
        cleaned_barcodes = {sample: "".join(re.findall(r"[A-Za-z]", barcode)) for sample, barcode in samples.items()}

        # Build a set of valid barcode sequences
        barcode_set = set(cleaned_barcodes.values())

        all_barcodes_including_hamming_distance = {}

        # Generate nearby sequences for each barcode
        for sample, cleaned_barcode in cleaned_barcodes.items():
            matches = set(gen_nearby_seqs(cleaned_barcode, barcode_set, max_hamming_distance))
            all_barcodes_including_hamming_distance[sample] = matches

        result[length] = all_barcodes_including_hamming_distance

    return result


def find_sample_for_read_index(index_str, sample_barcode_dict: dict) -> str:
    """
    Finds the sample name corresponding to a given barcode string.

    This function searches through a nested barcode dictionary to identify which sample
    the extracted barcode string belongs to.

    Args:
        index_str (str): The barcode string to search for.
        sample_barcode_dict (dict): A nested dictionary where outer keys are read lengths and inner
                             keys are sample names mapping to sets of barcode strings.

    Returns:
        str | undetermined: The sample name if a match is found; otherwise, undetermined.
    """
    if index_str is None or sample_barcode_dict is None:
        raise ValueError("Input is None.")

    if not isinstance(index_str, str):
        raise ValueError(f"{index_str} must be a string.")

    if not isinstance(sample_barcode_dict, dict):
        raise ValueError(f"{sample_barcode_dict} must be a dictionary.")

    for samples in sample_barcode_dict.values():
        for sample_name, barcodes in samples.items():
            if index_str in barcodes:
                return sample_name
    # If no match is found, return "undetermined"
    return "undetermined"


def demultiplex_fastq_by_barcode(fastq_file: str, samples_barcode_from_dict: dict, max_hamming_distance: int = 0, output_dir: str = ".") -> None:
    """
    Demultiplexes a FASTQ file by matching read indexes to known sample barcodes.

    This function extracts indexes from each read in a FASTQ file, compares them to a set
    of known sample barcodes (allowing for mismatches up to a given Hamming distance),
    and writes each read's sequence to a file named after the corresponding sample.
    Unmatched reads are written to an 'undetermined.txt' file.

    Args:
        fastq_file (str): Path to the FASTQ file containing sequencing reads.
        samples_barcode_from_csv (dict): Dictionary mapping sample names to sets of barcode sequences.
        max_hamming_distance (int, optional): Maximum number of mismatches allowed when comparing indexes.
                                              Defaults to 0 (exact match only).
        output_dir (str, optional): Directory where output files will be saved. Defaults to the current directory.

    Returns:
        None
    """

    # Group samples by index length, then generate sequences based on a pre-determined Hamming sequence
    grouped_samples_by_length = group_samples_by_index_length(samples_barcode_from_dict)
    all_barcodes_including_hamming_distance = generate_nearby_barcodes_by_length(grouped_samples_by_length, max_hamming_distance)

    file_handles = defaultdict(lambda: open(os.path.join(output_dir, "undetermined.txt"), "a", encoding="utf-8"))
    for sample in samples_barcode_from_dict:
        file_handles[sample] = open(os.path.join(output_dir, f"{sample}.txt"), "a", encoding="utf-8")

    # Read the FASTQ file and extract indexes, save all indexes in a list
    fastq = FastqFile(fastq_file)
    for name, seq, _ in fastq.open_read_iterator(as_string=True):
        # Extract the index from the read name and remove any non-alphabetic characters
        index = extract_index_from_header_illumina(name)
        index = re.sub(r"[^A-Za-z]", "", index)

        # Check if the index has any N's in it
        read_index_seq_list = []
        index_n_count = index.count('N')
        if index_n_count >= 1:
            barcode_set = {index}
            read_index_seq_list = gen_nearby_seqs(index, barcode_set, index_n_count)

        if read_index_seq_list is None:
            read_index_seq_list = [index]

        for seq in read_index_seq_list:
            # Check if the index is in the matches to any of the samples
            match_sample = find_sample_for_read_index(index, all_barcodes_including_hamming_distance)

            file_handles[match_sample].write(seq + "\n")

    for f in file_handles.values():
        f.close()


#############################################################################
# old approach
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
