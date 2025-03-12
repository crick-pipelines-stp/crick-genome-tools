"""
Class for generating vector core AAV report.
"""

import logging
from crick_genome_tools.reporting.reports.crick_report import CrickReport
from crick_genome_tools.reporting.tqc.plotly_charts import read_count_histogram, read_length_scatterplot, mqc_samtools_bar_plot, mqc_samtools_contig_bar_plot

import streamlit as st

log = logging.getLogger(__name__)

class VectorCoreAavReport(CrickReport):
    """
    Class for generating vector core AAV report.
    """

    def __init__(self, data_path, seq_mode):
        self.seq_mode = seq_mode
        super().__init__(data_path, "Vectorcore AAV Report")

    def generate_report(self, section_headers):
        section_headers = [
            "Read QC",
            "Contaminant Removal",
        ]
        super().generate_report(section_headers)

        # Activate current section
        self.activate_section()

    def activate_section(self):
        # Display the selected section content
        st.title(st.session_state.selected_section)

        if st.session_state.selected_section == "Read QC":
            self._read_qc_section()
        elif st.session_state.selected_section == "Contaminant Removal":
            self._contaminant_removal_section()

    def _read_qc_section(self):
        # Init
        dp = st.session_state.data_parser
        results_dict = dp.result_dict
        dataframe_dict = dp.dataframe_dict
        # st.write("This section shows read quality reporting.")

        # Create dropdown for selecting dataset
        selected_dataset = st.selectbox("Choose a sample:", list(results_dict.keys()))

        # Place charts and tables
        read_count_histogram(results_dict[selected_dataset]["toulligqc"])
        read_length_scatterplot(dataframe_dict[selected_dataset]["toulligqc"])

    def _contaminant_removal_section(self):
        # Init
        dp = st.session_state.data_parser

        # Prepare data
        host_df = dp.merged_dataframe_dict["samtools_host"]
        contam_df = dp.merged_dataframe_dict["samtools_contam"]
        contam_columns = contam_df.columns[1:].tolist()
        contam_columns.remove("Unmapped")
        host_df = host_df.reset_index(drop=True)
        combined_df = contam_df.reset_index(drop=True)
        combined_df.insert(2, "Host", host_df["Mapped"])
        combined_df.insert(1, "Total", contam_df.iloc[:, 1:].sum(axis=1))
        combined_columns = combined_df.columns[2:].tolist()
        combined_columns.remove("Unmapped")

        # Place charts and tables
        mqc_samtools_contig_bar_plot(combined_df, "Summary", combined_columns)
        mqc_samtools_bar_plot(dp.merged_dataframe_dict["samtools_host"], "Host Alignment")
        mqc_samtools_contig_bar_plot(contam_df, "Contaminent Alignment", contam_columns)
