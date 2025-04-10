"""
Manage SnapGene files.
"""

import csv


def convert_to_gff(input_path, output_path, contig_name, source="SnapGene"):
    """
    Convert a SnapGene file to GFF format.
    """
    with open(input_path, newline="", encoding="UTF-8") as infile, open(output_path, "w", newline="", encoding="UTF-8") as outfile:
        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t")

        for row in reader:
            if len(row) < 5:
                continue  # skip malformed lines

            attributes = "Note=" + row[0].strip()
            start_end = row[1].strip().split("..")
            start = start_end[0]
            end = start_end[1]
            strand_symbol = row[3].strip()
            feature_type = row[4].strip()

            strand = {"=>": "+", "==": "+", "<=": "-"}.get(strand_symbol, ".")

            gff_row = [
                contig_name,  # seqid
                source,  # source
                feature_type,  # type
                start,  # start
                end,  # end
                ".",  # score
                strand,  # strand
                ".",  # phase
                attributes,  # attributes
            ]
            writer.writerow(gff_row)
