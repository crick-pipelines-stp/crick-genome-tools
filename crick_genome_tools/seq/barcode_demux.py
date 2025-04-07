from crick_genome_tools.io.fastq_file import FastqFile


def extract_index_from_header(name: str) -> str:
    """
    Extract the index from the FASTQ read header.

    Args:
        name (str): The read name from the FASTQ file.

    Returns:
        str: The extracted index.
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
