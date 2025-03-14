"""
Class to manage reading of samtools files
"""


def parse_flagstat(flagstat_file):
    """
    Parse the flagstat file to extract total reads, mapped reads, and alignment rate.
    """
    with open(flagstat_file, "r", encoding="UTF-8") as f:
        lines = f.readlines()

    total_reads_line = next((line for line in lines if "total" in line.lower()), None)
    mapped_reads_line = next((line for line in lines if "mapped" in line.lower()), None)

    if total_reads_line is None or mapped_reads_line is None:
        raise ValueError("Required lines not found in the flagstat file.")

    # Extract total reads
    total_reads = int(total_reads_line.split("+")[0].strip())

    # Extract mapped reads
    mapped_reads = int(mapped_reads_line.split("+")[0].strip())

    # Calculate alignment rate from total and mapped reads
    if total_reads == 0:
        alignment_rate = 0
    else:
        alignment_rate = round(mapped_reads / total_reads * 100, 2)

    return total_reads, mapped_reads, alignment_rate
