"""
Utility functions for the blast application.
"""

import os
from typing import Dict
from collections import defaultdict

from crick_genome_tools.io.fasta import Fasta


def build_fasta_from_top_hits(fasta_path, blast_path) -> Dict[str, str]:
    """
    Build a fasta data structure from the top hits of a blast file.

    Args:
        fasta_path (str): Path to the reference fasta file.
        blast_path (str): Path to the blast file.

    Returns:
        Dict[str, str]: A dictionary where keys are the sequence identifiers and values are strings of the sequences.
    """
    # check if the fasta file exists
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"File does not exist: {fasta_path}")

    # check if the blast file exists
    if not os.path.exists(blast_path):
        raise FileNotFoundError(f"File does not exist: {blast_path}")

    # parse the blast file
    with open(blast_path, "r", encoding="UTF-8") as in_f:
        blast_lines = in_f.readlines()

    # if the blast file is empty, error out
    if len(blast_lines) == 0:
        raise ValueError("Error: No top hits detected in the blast file.")

    # parse the blast lines and group by query
    blast_hits = defaultdict(list)
    for line in blast_lines:
        parts = line.split("\t")
        query = parts[0]
        hit = parts[1]
        e_value = float(parts[10])
        coverage = float(parts[2])
        blast_hits[query].append((hit, e_value, coverage))

    # select the best hit for each query
    top_blast_hits = []
    for query, hits in blast_hits.items():
        best_hit = min(hits, key=lambda x: (x[1], -x[2]))  # lowest e-value, highest coverage
        top_blast_hits.append(best_hit[0])

    # parse the reference fasta file
    fasta_data = Fasta.read_fasta_file(fasta_path)

    # build the fasta data from the top hits
    fasta_top_hits = {}
    for hit in top_blast_hits:
        if hit in fasta_data:
            fasta_top_hits[hit] = fasta_data[hit]
        else:
            raise ValueError(f"Error: Hit {hit} not found in the reference fasta file.")

    return fasta_top_hits
