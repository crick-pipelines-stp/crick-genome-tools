"""
Contains access methods for using MutiQC to parse data.
"""

import logging
import multiqc
import pandas as pd


log = logging.getLogger(__name__)

def parse_samtools_stats(data_folder, sample_clean):
    """
    Parse samtools flagstats data.
    """
    # Parse data
    multiqc.parse_logs(data_folder, extra_fn_clean_exts=sample_clean)
    data = multiqc.get_general_stats_data()

    # Extract summary table
    rows = []
    for sample, metrics in data.items():
        mapped = metrics['Samtools: flagstat.mapped_passed']
        total = metrics['Samtools: flagstat.flagstat_total']
        unmapped = total - mapped
        rows.append({'Sample': sample, 'Total': total, 'Mapped': mapped, 'Unmapped': unmapped})

    df = pd.DataFrame(rows)
    df = df.sort_values(by=['Sample'])
    return df
