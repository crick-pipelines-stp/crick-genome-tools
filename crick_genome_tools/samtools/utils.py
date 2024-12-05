"""
Samtools utility functions.
"""

from collections import defaultdict


def count_table_from_pileup(pileup_path: str, output_path: str):
    """
    Generate a count table from a samtools mpileup file.

    Args:
        pileup_path (str): The path to the samtools mpileup file.
        output_path (str): The path to write the count table to.
    """
    # Read the pileup file
    with open(pileup_path, "r", encoding="UTF-8") as pileup_file:
        pileup_lines = pileup_file.readlines()

    # Parse the pilup lines
    data = {}
    for line in pileup_lines:
        fields = line.strip().split("\t")

        # Extract the fields
        contig = fields[0]
        position = int(fields[1])
        ref_base = fields[2]
        coverage = int(fields[3])
        bases = fields[4]
        qualities = fields[5]

        # Decode quality scores and get average
        qualities = [ord(quality) - 33 for quality in qualities]
        avg_quality = int(sum(qualities) / len(qualities))

        # Loop each base and count
        forward_count = 0
        rev_count = 0
        count_data = defaultdict(int)
        count_data["A"] = 0
        count_data["C"] = 0
        count_data["G"] = 0
        count_data["T"] = 0
        for base_char in bases:
            base = ref_base
            if base_char == ".":
                forward_count += 1
            elif base_char == ",":
                rev_count += 1
            elif base_char in "ACGT":
                base = base_char
            else:
                continue
            count_data[base] += 1

        # Turn count data to percentage rounded to 2 decimal places
        for base in count_data:
            count_data[base] = count_data[base] / coverage
            count_data[base] = round(count_data[base] * 100, 2)

        # Create data row
        data[position] = {
            "contig": contig,
            "ref_base": ref_base,
            "coverage": coverage,
            "avg_quality": avg_quality,
            "forward_count": forward_count,
            "rev_count": rev_count,
            "A": count_data["A"],
            "C": count_data["C"],
            "G": count_data["G"],
            "T": count_data["T"],
        }

    # Write the data to the output file
    with open(output_path, "w", encoding="UTF-8") as output_file:
        output_file.write("position\tcontig\tref\tcoverage\tavg_qual\tfwd_cnt\trev_cnt\tA%\tC%\tG%\tT%\n")
        for position, data in data.items():
            output_file.write(
                f"{position}\t{data['contig']}\t{data['ref_base']}\t{data['coverage']}\t{data['avg_quality']}\t{data['forward_count']}\t{data['rev_count']}\t{data['A']}\t{data['C']}\t{data['G']}\t{data['T']}\n"
            )
