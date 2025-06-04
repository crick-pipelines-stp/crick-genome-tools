"""
Contains access methods custom parsing of report data.
"""

import logging
import os

import pandas as pd

from crick_genome_tools.io.gzip_file import GzipFile


log = logging.getLogger(__name__)

# def to_cumulative_depth_curve(df, max_depth=100, min_percentage=0.02):
#     """
#     Convert to cumulative ≥depth percentages.
#     Ensures percentage values are normalized to sum to 1.0 before cumsum.
#     """
#     df = df.sort_values("Depth").copy()

#     total = df["Percentage"].sum()
#     if total == 0:
#         return df.head(0)  # return empty if no signal

#     df["Percentage"] = df["Percentage"] / total
#     df["Percentage"] = df["Percentage"][::-1].cumsum()[::-1]

#     df = df[df["Depth"] <= max_depth]
#     df = df[df["Percentage"] >= min_percentage]
#     return df

# def parse_mosdepth_globaldist(data_folder):
#     """
#     Parse mosdepth global dist data

#     contig	depth	percentage
#     P220_AAV_2100bp	185138	0.00
#     P220_AAV_2100bp	185130	0.00
#     P220_AAV_2100bp	185096	0.00
#     P220_AAV_2100bp	185095	0.00
#     P220_AAV_2100bp	185087	0.00
#     P220_AAV_2100bp	185034	0.00
#     P220_AAV_2100bp	185029	0.00
#     P220_AAV_2100bp	185027	0.00
#     P220_AAV_2100bp	185016	0.00
#     P220_AAV_2100bp	185013	0.00
#     P220_AAV_2100bp	185012	0.01

#     Parses all *.mosdepth.global.dist.txt files in a folder into a single long-form dataframe
#     with columns: Sample, Contig, Depth, Percentage
#     """
#     all_dfs = []

#     for file in os.listdir(data_folder):
#         if not file.endswith(".mosdepth.global.dist.txt"):
#             continue

#         sample_id = file.replace(".mosdepth.global.dist.txt", "")
#         file_path = os.path.join(data_folder, file)

#         with open(file_path, "r", encoding="UTF-8") as f:
#             rows = []
#             for line in f:
#                 contig, depth, percentage = line.strip().split("\t")
#                 rows.append((sample_id, contig, int(depth), float(percentage)))

#         df = pd.DataFrame(rows, columns=["Sample", "Contig", "Depth", "Percentage"])
#         all_dfs.append(df)

#     combined_df = pd.concat(all_dfs, ignore_index=True)
#     df_cum = combined_df.groupby(["Sample", "Contig"], group_keys=False).apply(
#         lambda d: to_cumulative_depth_curve(d, max_depth=100)
#     )
#     return df_cum


def parse_mosdepth_per_base(data_folder):
    """
    Parse mosdepth per base data

    P220_AAV_2100bp	0	1	129208
    P220_AAV_2100bp	1	2	130866
    P220_AAV_2100bp	2	3	132201
    P220_AAV_2100bp	3	4	133771
    P220_AAV_2100bp	4	5	136160
    """
    data = {}

    # Iterate over all files in the folder
    for file in os.listdir(data_folder):
        if not file.endswith(".per-base.bed.gz"):
            continue

        # Add sample id
        sample_id = file.replace(".per-base.bed.gz", "")
        if sample_id not in data:
            data[sample_id] = {}

        # Read gzipped bed file
        bed_file = GzipFile(os.path.join(data_folder, file))
        per_contig = {}
        for line in bed_file.open_read_iterator(as_string=True):
            fields = line.strip().split("\t")
            contig = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            depth = int(fields[3])

            if contig not in per_contig:
                per_contig[contig] = []

            per_contig[contig].extend((pos, depth) for pos in range(start, end))

        # For each contig, convert to DataFrame and fill gaps
        for contig, entries in per_contig.items():
            df = pd.DataFrame(entries, columns=["Position", "Depth"])
            min_pos = df["Position"].min()
            max_pos = df["Position"].max()
            full_range = pd.DataFrame({"Position": range(min_pos, max_pos + 1)})
            df = full_range.merge(df, on="Position", how="left").fillna(0)
            df["Depth"] = df["Depth"].astype(int)
            data[sample_id][contig] = df

    # Sort dict by sample id
    data = dict(sorted(data.items()))

    return data
