import os
import re
from collections import defaultdict
import numpy as np

from crick_genome_tools.io.fastq_file import FastqFile


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

def vectorised_find_closest_match(barcode_dict: dict[str, str], query_seq: str, max_hamming: int) -> str:
    """
    Finds the sample name with the closest matching barcode to the query sequence.
    If sequences are of unequal length, comparison is done up to the length of the shorter one.

    Args:
        barcode_dict (dict): Dictionary mapping sample names to barcode strings.
        query_seq (str): Sequence to compare against.
        max_hamming (int): Maximum allowed Hamming distance.

    Returns:
        str: The name of the best matching sample, or "undetermined" if no match within max_hamming.
    """
    best_sample = "undetermined"
    min_distance = float("inf")

    query_len = len(query_seq)
    query_array = np.array(list(query_seq), dtype="U1")

    for sample, barcode in barcode_dict.items():
        compare_len = min(query_len, len(barcode))
        if compare_len == 0:
            continue  # skip empty sequences

        barcode_array = np.array(list(barcode[:compare_len]), dtype="U1")
        trimmed_query = query_array[:compare_len]

        mismatches = barcode_array != trimmed_query
        distance = np.sum(mismatches)

        if distance < min_distance and distance <= max_hamming:
            min_distance = distance
            best_sample = sample
            if distance == 0:
                break  # perfect match

    return best_sample


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
    for group in grouped_samples_by_length.values():
        for sample in group:
            group[sample] = re.sub(r"[^A-Za-z]", "", group[sample])

    file_handles = defaultdict(lambda: open(os.path.join(output_dir, "undetermined.txt"), "a", encoding="utf-8"))
    for sample in samples_barcode_from_dict:
        file_handles[sample] = open(os.path.join(output_dir, f"{sample}.txt"), "a", encoding="utf-8")

    # Read the FASTQ file and extract indexes, save all indexes in a list
    fastq = FastqFile(fastq_file)
    for name, seq, _ in fastq.open_read_iterator(as_string=True):
        # Extract the index from the read name and remove any non-alphabetic characters
        index = extract_index_from_header_illumina(name)
        index = re.sub(r"[^A-Za-z]", "", index)

        # Search for the closest match in the grouped samples from the longest to the shortest indexes
        match = "undetermined"
        for length in sorted(grouped_samples_by_length.keys(), reverse=True):
            trimmed_index = index[:length]
            group = grouped_samples_by_length[
                length
            ]  # this trims the index at the end of the sequence, ingores the possibility of the sequence being 2 indexes merged into one
            print(group)

            # Check if the read index matches to any of the samples
            match = vectorised_find_closest_match(group, trimmed_index, max_hamming_distance)
            print(match)
            if match != "undetermined":
                break

        # Write the sequence to the appropriate file
        file_handles[match].write(seq + "\n")

    for f in file_handles.values():
        f.close()
