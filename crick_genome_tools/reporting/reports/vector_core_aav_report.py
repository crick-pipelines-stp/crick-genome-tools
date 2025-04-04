"""
Class for generating vector core AAV report.
"""

# pylint disable=missing-function-docstring,missing-class-docstring

import logging
import json

import streamlit as st

from crick_genome_tools.reporting.reports.crick_report import CrickReport
from crick_genome_tools.reporting.tqc.plotly_charts import read_count_histogram, read_length_scatterplot, mqc_samtools_bar_plot, mqc_samtools_contig_bar_plot, coverage_plot
from crick_genome_tools.io.fasta import Fasta


log = logging.getLogger(__name__)

class VectorCoreAavReport(CrickReport):
    """
    Class for generating vector core AAV report.
    """

    def __init__(self, run_id, data_path = None, data_obj = None, tmp_dir = None, jbrowse_component = None, app_url = "localhost:8501"):
        super().__init__("Vectorcore AAV Report", data_path, data_obj, tmp_dir)
        self.run_id = run_id
        self.jbrowse_component = jbrowse_component
        self.app_url = app_url

    def generate_report(self, section_headers = []):
        section_headers = [
            "Pipeline Summary",
            "Read QC",
            "Contaminant Removal",
            "Alignment",
            "Coverage",
            "Genome Viewer",
        ]
        super().generate_report(section_headers)
        # st.subheader(self.run_id)

        # Activate current section
        self.activate_section()

    def activate_section(self):
        # Display the selected section content
        st.title(st.session_state.selected_section)
        dp = st.session_state.data_parser

        if st.session_state.selected_section == "Pipeline Summary":
            self.summary_section(dp)
        elif st.session_state.selected_section == "Read QC":
            self.read_qc_section(dp)
        elif st.session_state.selected_section == "Contaminant Removal":
            self.contaminant_removal_section(dp)
        elif st.session_state.selected_section == "Alignment":
            self.alignment_section(dp)
        elif st.session_state.selected_section == "Coverage":
            self.coverage_section(dp)
        elif st.session_state.selected_section == "Genome Viewer":
            self.genome_viewer_section(dp)

    def render_table(self, param_dict):
        lines = ["| Parameter | Value |", "|---|---|"]
        for key, value in param_dict.items():
            lines.append(f"| {key} | {value} |")
        return "\n".join(lines)

    def summary_section(self, dp):
        for section, section_params in dp.summary_data.items():
            with st.expander(section, expanded=True):
                st.markdown(self.render_table(section_params), unsafe_allow_html=True)

    def read_qc_section(self, dp):
        # Init
        results_dict = dp.result_dict
        dataframe_dict = dp.dataframe_dict
        # st.write("This section shows read quality reporting.")

        # Create dropdown for selecting dataset
        selected_dataset = st.selectbox("Choose a sample:", list(results_dict.keys()))

        # Place charts and tables
        read_count_histogram(results_dict[selected_dataset]["toulligqc"])
        read_length_scatterplot(dataframe_dict[selected_dataset]["toulligqc"])

    def contaminant_removal_section(self, dp):
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

    def alignment_section(self, dp):
        # Place charts and tables
        mqc_samtools_bar_plot(dp.merged_dataframe_dict["samtools_align"], "AAV Alignment")

    def coverage_section(self, dp):
        # Get data
        coverage_data = dp.dataframe_dict["coverage_per_base"]

        # Create dropdown for selecting dataset
        selected_dataset = st.selectbox("Choose a sample:", list(coverage_data.keys()))

        # Place chart for each contig
        for contig, df in coverage_data[selected_dataset].items():
            coverage_plot(df, contig)

    def genome_viewer_section(self, dp):
        # Create dropdown for selecting dataset
        selected_dataset = st.selectbox("Choose a sample:", list(dp.result_dict.keys()))

        # Write data to static folder
        ref_path = self.tmp_dir + "_" + selected_dataset + ".fasta"
        with open(ref_path, "wb") as f:
            for line in dp.result_dict[selected_dataset]["ref"]:
                f.write(line.encode("utf-8"))
        index_path = self.tmp_dir + "_" + selected_dataset + ".fai"
        with open(index_path, "wb") as f:
            for line in dp.result_dict[selected_dataset]["fai"]:
                f.write(line.encode("utf-8"))

        # Construct Uris
        base_uri = self.app_url + "/app/static/tmp/" + self.tmp_dir.split("/")[-1] + "_" + selected_dataset
        fasta_uri = base_uri + ".fasta"
        fai_uri = base_uri + ".fai"

        # Read the fasta file
        fasta_data = Fasta.read_fasta_file(ref_path)
        contigs = list(fasta_data.keys())

        # fasta_data = fasta_data[list(fasta_data.keys())[0]]
        # contig_length = len(fasta_data)

        # Build the jbrowse config
        jbrowse_config = {
            "assembly": {
                "name": f"{contigs[0]}",
                "sequence": {
                    "type": "ReferenceSequenceTrack",
                    "trackId": f"{contigs[0]}_refseq",
                    "adapter": {
                        "type": "IndexedFastaAdapter",
                        "fastaLocation": {
                            "uri": f"{fasta_uri}",
                            "locationType": "UriLocation"
                        },
                        "faiLocation": {
                            "uri":f"{fai_uri}",
                            "locationType": "UriLocation"
                        },
                    }
                }
            },
            "tracks": [],
            "defaultSession": {
                "name": "Default session",
                "margin": 0,
                "view": {
                    "id": "linearGenomeView",
                    "type": "LinearGenomeView",
                    "init": {
                        "assembly": 'hg38',
                        "loc": '1..100',
                        "tracks": []
                    }
                },
            }
        }
        # print(json.dumps(jbrowse_config, indent=4))

        if self.jbrowse_component is not None:
            self.jbrowse_component("aav_viewer", config=jbrowse_config, height=800)
