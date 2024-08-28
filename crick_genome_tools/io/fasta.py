"""
Class to manage reading of a FASTA file
"""

import logging
import os
from typing import Dict, List


log = logging.getLogger(__name__)


class Fasta:
    """
    Class to manage reading of a FASTA file to a data structure
    """

    @staticmethod
    def read_fasta_file(path: str) -> Dict[str, str]:
        """
        Read the FASTA file from path into a dictionary of strings

        Args:
            path (str): Path to FASTA file.

        Returns:
            Dict[str, str]: A dictionary where keys are the sequence identifiers and values are strings of the sequences.
        """

        # Check for path not exists
        if not os.path.exists(path):
            raise FileNotFoundError(f"File does not exist: {path}")

        # Init
        sequences: Dict[str, str] = {}
        current_tag = None
        current_seq: List[str] = []

        with open(path, "r", encoding="UTF-8") as file:
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    if current_tag is not None:
                        sequences[current_tag] = "".join(current_seq)
                    # Start a new sequence
                    current_tag = line[1:].strip()  # Skip the '>' and take the rest as the identifier
                    current_seq = []
                else:
                    # Append this line of sequence
                    current_seq.extend(line.upper())

            # Add final seq to dict
            if current_tag is not None:
                sequences[current_tag] = "".join(current_seq)

        # Check if the dictionary is empty and raise an error if it is
        if not sequences:
            raise ValueError(f"The FASTA file is empty or incorrectly formatted: {path}")

        return sequences

    @staticmethod
    def read_fasta_directory(path: str) -> Dict[str, str]:
        """
        Read a directory of fasta files

        Args:
            path (str): Path to FASTA directory.

        Returns:
            Dict[str, str]: A dictionary where keys are the sequence identifiers and values are strings of the sequences.
        """

        # Init
        sequences: Dict[str, str] = {}

        # Check exists
        if not os.path.isdir(path):
            raise FileNotFoundError("Directory does not exist")

        # Read files
        files = [
            os.path.join(path, file)
            for file in os.listdir(path)
            if os.path.isfile(os.path.join(path, file)) and (file.endswith(".fasta") or file.endswith(".fa"))
        ]
        for file in files:
            sequences.update(Fasta.read_fasta_file(file))

        return sequences

    @staticmethod
    def write_fasta_file(fasta_data: Dict[str, str], path: str, column_width: int = 70):
        """
        Write a dictionary of fasta data to a file

        Args:
            fasta_data (Dict[str, str]): Dictionary of fasta data.
            path (str): Path to write the fasta file to.
        """

        with open(path, "w", encoding="UTF-8") as file:
            for tag, seq in fasta_data.items():
                file.write(f">{tag}\n")
                for i in range(0, len(seq), column_width):
                    file.write(seq[i:i+column_width] + "\n")
