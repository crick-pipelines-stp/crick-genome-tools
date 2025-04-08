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


def read_fastq(self, fastq_file):
    """
    Read a FASTQ file and return a list of sequences.

    Args:
        file_path (str): Path to the FASTQ file.

    Returns:
        list: List of sequences.
    """
    fastq = FastqFile(fastq_file)

    for name, seq, qual in fastq.open_read_iterator(as_string=True):
        self.extract_index_from_header(name)

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
