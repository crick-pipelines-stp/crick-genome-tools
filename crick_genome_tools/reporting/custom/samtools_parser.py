"""
Contains access methods custom parsing of report data.
"""

import os
import logging


log = logging.getLogger(__name__)

def parse_samtools_flagstat(data_folder, sample_clean):
    """
    Parse samtools flagstats data.

    123018 + 0 in total (QC-passed reads + QC-failed reads)
    122215 + 0 primary
    180 + 0 secondary
    623 + 0 supplementary
    0 + 0 duplicates
    0 + 0 primary duplicates
    1482 + 0 mapped (1.20% : N/A)
    679 + 0 primary mapped (0.56% : N/A)
    0 + 0 paired in sequencing
    0 + 0 read1
    0 + 0 read2
    0 + 0 properly paired (N/A : N/A)
    0 + 0 with itself and mate mapped
    0 + 0 singletons (N/A : N/A)
    0 + 0 with mate mapped to a different chr
    0 + 0 with mate mapped to a different chr (mapQ>=5)
    """
    # Init
    data = {}

    # Find files in folder with .flagstat extension
    file_list = []
    for file in os.listdir(data_folder):
        if file.endswith(".flagstat"):
            file_list.append(file)

    # Cycle through files
    for file_name in file_list:
        # Clean file name
        sample_name = file_name.replace(sample_clean, "").replace(".flagstat", "")
        data[sample_name] = {}
        log.debug(f"Processing {file_name} for {sample_name}.")

        # Parse data
        with open(os.path.join(data_folder, file_name), "r", encoding="UTF-8") as file:
            for line in file:
                if "total (QC-passed reads + QC-failed reads)" in line:
                    data[sample_name]["total_aligned"] = int(line.split()[0])
                elif "primary mapped" in line:
                    data[sample_name]["primary_mapped"] = int(line.split()[0])
                    break
                elif "primary duplicates" in line:
                    data[sample_name]["primary_duplicates"] = int(line.split()[0])
                elif "primary" in line:
                    data[sample_name]["primary"] = int(line.split()[0])
                elif "secondary" in line:
                    data[sample_name]["secondary"] = int(line.split()[0])
                elif "supplementary" in line:
                    data[sample_name]["supplementary"] = int(line.split()[0])
                elif "duplicates" in line:
                    data[sample_name]["duplicates"] = int(line.split()[0])
                elif "mapped" in line:
                    data[sample_name]["mapped"] = int(line.split()[0])
    return data
