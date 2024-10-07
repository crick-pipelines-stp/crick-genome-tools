"""
Class to manage reading of samtools files
"""

def parse_flagstat(flagstat_file):
    """
    Parse the flagstat file to extract total reads, mapped reads, and alignment rate.
    """
    with open(flagstat_file, 'r', encoding="UTF-8") as f:
        lines = f.readlines()

    total_reads_line = lines[0]
    mapped_reads_line = lines[4]

    # Extract total reads
    total_reads = int(total_reads_line.split('+')[0].strip())
    # Extract mapped reads
    mapped_reads = int(mapped_reads_line.split('+')[0].strip())
    # Extract alignment rate
    alignment_rate_str = mapped_reads_line.strip().split('(')[1].split('%')[0]
    alignment_rate = float(alignment_rate_str)

    return total_reads, mapped_reads, alignment_rate
