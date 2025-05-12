import os
import re
from collections import defaultdict

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


import numpy as np

def vectorised_hamming(query_seq: str, barcodes: list[str], max_hamming: int):
    """
    Finds the barcode with the smallest Hamming distance using NumPy vectorization.

    Args:
        query_seq (str): The sequence to compare.
        barcodes (list[str]): List of barcode strings of the same length as query_seq.
        max_hamming (int): Maximum allowed Hamming distance.

    Returns:
        tuple: (best matching barcode, hamming distance), or (None, None) if no match.
    """
    if not barcodes:
        return None, None

    # Convert query and barcodes into numpy arrays of shape (num_barcodes, seq_length)
    barcode_array = np.array([list(b) for b in barcodes], dtype='U1')  # shape (n, k)
    query_array = np.array(list(query_seq), dtype='U1')                # shape (k,)

    # Broadcast and compare: this returns a boolean array (n, k)
    mismatches = barcode_array != query_array

    # Sum mismatches (hamming distance) across each barcode (axis=1)
    distances = np.sum(mismatches, axis=1)

    # Find the smallest distance under the threshold
    min_dist = np.min(distances)
    if min_dist <= max_hamming:
        best_idx = np.argmin(distances)
        return barcodes[best_idx], min_dist
    else:
        return None, None



def find_closest_match(barcode_dict: dict, seq: str, max_hamming: int) -> str:
    """
    Finds the sample and barcode with the smallest Hamming distance to the given sequence.

    Args:
        barcode_dict (dict): Dictionary mapping sample names to barcode strings.
        seq (str): The sequence to compare against the barcodes.
        max_hamming (int): Maximum allowed Hamming distance.

    Returns:
        str or None: Returns sample_name with the smallest
            Hamming distance, or None if no barcode is within max_hamming.
    """
    if not isinstance(seq, str):
        raise ValueError(f"{seq} must be a string.")

    if not isinstance(barcode_dict, dict):
        raise ValueError(f"{barcode_dict} must be a dictionary.")

    if not isinstance(max_hamming, int):
        raise ValueError(f"{max_hamming} must be an integer.")

    if max_hamming < 0:
        raise ValueError(f"{max_hamming} must be a non-negative integer.")

    best_sample = "undetermined"
    min_distance = float("inf")

    for sample, barcodes in barcode_dict.items():
        match, dist = vectorised_hamming(seq, barcodes, max_hamming)
        if match and dist < min_distance:
            best_sample = sample
            min_distance = dist
            if dist == 0:  # perfect match
                break

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


import numpy as np

def vectorized_find_closest_sample(barcode_dict: dict[str, str], query_seq: str, max_hamming: int) -> str:
    """
    Finds the sample name with the closest matching barcode to the query sequence using NumPy vectorization.
    
    Args:
        barcode_dict (dict): Dictionary mapping sample names to a single barcode string each.
        query_seq (str): Sequence to compare against.
        max_hamming (int): Maximum allowed Hamming distance.
    
    Returns:
        str: The name of the best matching sample, or "undetermined" if no match within max_hamming.
    """
    best_sample = "undetermined"
    min_distance = float("inf")

    # Filter only barcodes of matching length
    matching = {k: v for k, v in barcode_dict.items() if len(v) == len(query_seq)}
    if not matching:
        return best_sample

    # Convert barcodes to a 2D NumPy array
    sample_names = list(matching.keys())
    barcodes = list(matching.values())
    barcode_array = np.array([list(b) for b in barcodes], dtype="U1")
    query_array = np.array(list(query_seq), dtype="U1")

    # Compute Hamming distances
    mismatches = barcode_array != query_array
    distances = np.sum(mismatches, axis=1)

    # Find best match
    min_idx = np.argmin(distances)
    min_dist = distances[min_idx]

    if min_dist <= max_hamming:
        best_sample = sample_names[min_idx]

    return best_sample



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
            match = vectorized_find_closest_sample(group, trimmed_index, max_hamming_distance)
            print(match)
            if match != "undetermined":
                break

        # Write the sequence to the appropriate file
        file_handles[match].write(seq + "\n")

    for f in file_handles.values():
        f.close()
