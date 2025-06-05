import gzip
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


def group_samples_by_index_length(sample_index_dict: dict) -> list:
    """
    Organizes sample barcode data by the lengths of their first and last components.

    This function processes a dictionary of sample barcode entries, where each value
    can be either a string (representing the barcode) or a dictionary containing keys
    such as "barcode", "index", and optionally "index2". The barcode string is split
    into components using non-alphabetic characters as delimiters. The lengths of the
    first and last components are used to create a tuple key, which maps to a nested
    dictionary of sample names and their corresponding barcode strings.

    Args:
        sample_index_dict (dict): A dictionary where keys are sample names and values
            are either:
                - A string representing a barcode.
                - A dictionary that may include "barcode", "index", and optionally "index2" keys.

    Returns:
        dict: A nested dictionary where each key is a tuple of the lengths of the first
            and last barcode components, and each value is a dictionary mapping sample
            names to their raw barcode strings.

    Raises:
        TypeError: If a sample value is neither a string nor a dictionary, or if a barcode
            extracted from a sample is not a string.
    """
    # Check if the input is a dictionary
    if not isinstance(sample_index_dict, dict):
        raise TypeError(f"{sample_index_dict} must be a dictionary.")

    result = defaultdict(dict)
    # Calculate barcode lengths and organize samples by these lengths within `result`
    for sample, value in sample_index_dict.items():
        if isinstance(value, str):
            # if sample has a string value, treat it as a barcode
            barcode = value

        elif isinstance(value, dict):
            # if sample has a dict value, look for "barcode" or "index" keys
            if "barcode" in value:
                barcode = value["barcode"]
            elif "index" in value:
                if "index2" in value:
                    # if both index and index2 are present, use both entries as the barcode
                    barcode = value["index"] + "," + value["index2"]
                else:
                    # if only index is present, use it as the barcode
                    barcode = value["index"]
            else:
                barcode = ""
        else:
            raise TypeError(f"Sample '{sample}' has an unsupported value type: {type(value).__name__}")

        if not isinstance(barcode, str):
            raise TypeError(f"Barcode for sample '{sample}' must be a string.")

        # Split the barcode into parts based on non-alphabetic characters
        parts = re.split(r"[^A-Za-z]+", barcode)

        # Determine the lengths of the first and last parts
        first_len = len(parts[0]) if parts else 0
        last_len = len(parts[-1]) if len(parts) > 1 else 0 if not parts else 0 if len(parts) == 1 else len(parts[-1])

        key = (first_len, last_len)
        result[key][sample] = barcode

    return dict(result)


def hamming_distance(seq1, seq2) -> int:
    """Computes the Hamming distance between two sequences."""
    if seq1 is None or seq2 is None:
        raise ValueError("Input sequences cannot be None.")

    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


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

    best_match = "undetermined"
    min_distance = float("inf")

    for sample, barcode in barcode_dict.items():
        dist = hamming_distance(seq, barcode)
        # If the strings are a perfect match, skip further comparisons
        if dist == 0:
            return sample

        if dist < min_distance and dist <= max_hamming:
            min_distance = dist
            best_match = sample

    return best_match


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

def index_to_match_key(index: str, grouped_samples_by_length: dict, max_hamming_distance: int) -> tuple:
    # Extract the index from the read name and remove any non-alphabetic characters
    index = extract_index_from_header_illumina(index)
    index = re.sub(r"[^A-Za-z]", " ", index)

    # Search for the closest match in the grouped samples from the longest to the shortest indexes
    match = "undetermined"
    for length_key in sorted(grouped_samples_by_length.keys(), key=custom_priority_by_length_sort_key):
        group = grouped_samples_by_length[length_key]

        # trimming differes depending on whether it's a single or dual index
        length = sum(length_key)
        trimmed_index = trim_merge_string(index, length)

        # Check if the read index matches to any of the samples
        match = find_closest_match(group, trimmed_index, max_hamming_distance)

        # Stop searching if a match with a defined sample is found
        if match != "undetermined":
            break
    
    return match, trimmed_index


def demultiplex_fastq_by_barcode(samples_barcode_from_dict: dict, fastq_file_R1: str, max_hamming_distance: int = 0, output_dir: str = ".", fastq_file_R2: str = None) -> dict:
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
    ## Group samples by index length, then generate sequences based on a pre-determined Hamming sequence
    grouped_samples_by_length = group_samples_by_index_length(samples_barcode_from_dict)
    for group in grouped_samples_by_length.values():
        for sample in group:
            # merge barcodes into an individual string for demux processing
            group[sample] = re.sub(r"[^A-Za-z]", "", group[sample])

    ## Compare the all the barcodes in each group against each other to find the Hamming distance for each pair compared
    grouped_sample_by_length_hamming_value = {}
    grouped_sample_by_length_hamming_value = {
        length: crosscheck_barcode_proximity(samples)
        for length, samples in grouped_samples_by_length.items()
        if len(samples) > 1  # Skip groups with only one entry
    }
    # find the minimum hamming distance for each index-length group
    min_hamming_distances_by_length = find_min_hamming_distances(grouped_sample_by_length_hamming_value)
    # check if the minimum hamming distance is above the threshold max hamming distance
    # returns ValueError if any group has a minimum Hamming distance less than max_hamming
    assert_min_hamming_above_threshold(min_hamming_distances_by_length, max_hamming_distance)

    ## Create a fastq file for each sample + an "undetermined" file for unassigned reads
    file_handles_R1 = {sample: gzip.open(os.path.join(output_dir, f"{sample}_R1.fastq.gz"), "ab") for sample in samples_barcode_from_dict}
    file_handles_R1["undetermined"] = gzip.open(os.path.join(output_dir, "undetermined_R1.fastq.gz"), "ab")

    if fastq_file_R2:
        file_handles_R2 = {sample: gzip.open(os.path.join(output_dir, f"{sample}_R2.fastq.gz"), "ab") for sample in samples_barcode_from_dict}
        file_handles_R2["undetermined"] = gzip.open(os.path.join(output_dir, "undetermined_R2.fastq.gz"), "ab")

    # Keep track of each read assigned to each sample
    sample_assigned_read = defaultdict(list)  # this one keeps track of sample+read index information
    sample_count = {}  # this one keeps track of how many reads are assigned to each sample

    ## Read the FASTQ file and extract indexes, save all indexes in a list
    fastq_1 = FastqFile(fastq_file_R1)
    if fastq_file_R2:
        fastq_2 = FastqFile(fastq_file_R2)
        fastq_2_iter = fastq_2.open_read_iterator(as_string=True)

    for name, seq, qual in fastq_1.open_read_iterator(as_string=True):

        match, trimmed_index = index_to_match_key(name, grouped_samples_by_length, max_hamming_distance)

        # Assign the read to the matched sample
        sample_assigned_read[match].append([name, trimmed_index])

        # Write the sequence to the appropriate file
        FastqFile.write_read(file_handles_R1[match], name, seq, qual)

        if fastq_file_R2:
            # Read the corresponding R2 read
            name_R2, seq_R2, qual_R2 = next(fastq_2_iter)

            # name_R2, seq_R2, qual_R2 = next(fastq_2.open_read_iterator(as_string=True))
            # Write the R2 read to the appropriate file
            FastqFile.write_read(file_handles_R2[match], name_R2, seq_R2, qual_R2)
            

        # # Extract the index from the read name and remove any non-alphabetic characters
        # index = extract_index_from_header_illumina(name)
        # index = re.sub(r"[^A-Za-z]", " ", index)

        # # Search for the closest match in the grouped samples from the longest to the shortest indexes
        # match = "undetermined"
        # for length_key in sorted(grouped_samples_by_length.keys(), key=custom_priority_by_length_sort_key):
        #     group = grouped_samples_by_length[length_key]

        #     # trimming differes depending on whether it's a single or dual index
        #     length = sum(length_key)
        #     trimmed_index = trim_merge_string(index, length)

        #     # Check if the read index matches to any of the samples
        #     match = find_closest_match(group, trimmed_index, max_hamming_distance)

        #     # Stop searching if a match with a defined sample is found
        #     if match != "undetermined":
        #         break

        # # Assign the read to the matched sample
        # sample_assigned_read[match].append([name, trimmed_index])

        # # Write the sequence to the appropriate file
        # FastqFile.write_read(file_handles[match], name, seq, qual)

    for sample, reads in sample_assigned_read.items():
        sample_count[sample] = len(reads)

    for f in file_handles_R1.values():
        f.close()

    if fastq_file_R2:
        for f in file_handles_R2.values():
            f.close()

    return sample_count
