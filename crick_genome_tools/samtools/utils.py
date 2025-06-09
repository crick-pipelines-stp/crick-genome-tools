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

    data = {}

    for line in pileup_lines:
        fields = line.strip().split("\t")
        if len(fields) < 6:
            continue

        contig = fields[0]
        position = int(fields[1])
        ref_base = fields[2].upper()
        coverage = int(fields[3])
        bases = fields[4]
        qualities = fields[5]

        if qualities:
            qual_vals = [ord(q) - 33 for q in qualities]
            avg_quality = int(sum(qual_vals) / len(qual_vals))
        else:
            avg_quality = 0

        forward_count = 0
        rev_count = 0
        count_data = defaultdict(int, {"A": 0, "C": 0, "G": 0, "T": 0})

        i = 0
        while i < len(bases):
            char = bases[i]

            if char == '^':
                i += 2  # skip ^ and mapping quality
                continue
            elif char == '$':
                i += 1
                continue
            elif char in '+-':
                i += 1
                indel_len_str = ''
                while i < len(bases) and bases[i].isdigit():
                    indel_len_str += bases[i]
                    i += 1
                indel_len = int(indel_len_str) if indel_len_str else 0
                i += indel_len
                continue
            elif char == '.':
                count_data[ref_base] += 1
                forward_count += 1
            elif char == ',':
                count_data[ref_base] += 1
                rev_count += 1
            elif char.upper() in "ACGT":
                count_data[char.upper()] += 1
                if char.isupper():
                    forward_count += 1
                else:
                    rev_count += 1
            # skip anything else
            i += 1

        if coverage > 0:
            for base in "ACGT":
                count_data[base] = round((count_data[base] / coverage) * 100, 2)
        else:
            for base in "ACGT":
                count_data[base] = 0.0

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

    with open(output_path, "w", encoding="UTF-8") as output_file:
        output_file.write("position\tcontig\tref\tcoverage\tavg_qual\tfwd_cnt\trev_cnt\tA%\tC%\tG%\tT%\n")
        for position in sorted(data.keys()):
            d = data[position]
            output_file.write(
                f"{position}\t{d['contig']}\t{d['ref_base']}\t{d['coverage']}\t{d['avg_quality']}\t{d['forward_count']}\t{d['rev_count']}\t{d['A']}\t{d['C']}\t{d['G']}\t{d['T']}\n"
            )
