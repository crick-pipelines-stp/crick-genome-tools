"""
Class for generating a generic viral genomics report.
"""

# pylint disable=missing-function-docstring,missing-class-docstring,missing-module-docstring,missing-function-docstring

import json
import logging

import pandas as pd
import streamlit as st
from pandas.api.types import is_categorical_dtype, is_datetime64_any_dtype, is_numeric_dtype

from crick_genome_tools.io.fasta import Fasta
from crick_genome_tools.reporting.reports.crick_report import CrickReport
from crick_genome_tools.reporting.tqc.plotly_charts import (
    coverage_plot,
    mqc_samtools_bar_plot,
    mqc_samtools_contig_bar_plot,
    read_count_histogram,
    read_length_scatterplot,
    truncation_barplot,
    truncation_scatterplot,
)


log = logging.getLogger(__name__)


class ViralGenomicsReport(CrickReport):
    """
    Base class for generating a generic viral genomics report.
    """

    def __init__(self, run_id, data_path=None, data_obj=None, tmp_dir=None, jbrowse_component=None, app_url="localhost:8501"):
        super().__init__("Viral Genomics Report", data_path, data_obj, tmp_dir)
        self.run_id = run_id
        self.jbrowse_component = jbrowse_component
        self.app_url = app_url

    def _scan_sections(self):
        """
        Scan the sections of the report and return a list of section headers.
        """

        dp = st.session_state.data_parser
        headers = []
        if dp.summary_data is not None:
            headers.append("Pipeline Summary")
        if "samplesheet" in dp.merged_dataframe_dict:
            headers.append("Samplesheet")
        if "toulligqc" in dp.result_dict:
            headers.append("Read QC")
        if "samtools_host" in dp.merged_dataframe_dict or "samtools_contam" in dp.merged_dataframe_dict:
            headers.append("Contaminant Removal")
        if "samtools_align" in dp.merged_dataframe_dict:
            headers.append("Alignment")
        if "coverage_per_base" in dp.dataframe_dict:
            headers.append("Coverage")
        if "consensus" in dp.result_dict:
            headers.append("Consensus")
        if "truncation" in dp.result_dict:
            headers.append("Truncation")
        if "reference" in dp.result_dict:
            headers.append("Genome Viewer")
        if "variants" in dp.dataframe_dict:
            headers.append("Variant Viewer")
        return headers

    def generate_report(self, section_headers=[]):
        section_headers = self._scan_sections()
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
        elif st.session_state.selected_section == "Samplesheet":
            self.samplesheet_section(dp)
        elif st.session_state.selected_section == "Read QC":
            self.nanopore_read_qc_section(dp)
        elif st.session_state.selected_section == "Contaminant Removal":
            self.contaminant_removal_section(dp)
        elif st.session_state.selected_section == "Alignment":
            self.alignment_section(dp)
        elif st.session_state.selected_section == "Coverage":
            self.coverage_section(dp)
        elif st.session_state.selected_section == "Truncation":
            self.truncation_section(dp)
        elif st.session_state.selected_section == "Consensus":
            self.consensus_section(dp)
        elif st.session_state.selected_section == "Genome Viewer":
            self.genome_viewer_section(dp)
        elif st.session_state.selected_section == "Variant Viewer":
            self.variant_viewer_section(dp)

    def render_table(self, param_dict):
        lines = ["| Parameter | Value |", "|---|---|"]
        for key, value in param_dict.items():
            lines.append(f"| {key} | {value} |")
        return "\n".join(lines)

    def summary_section(self, dp):
        for section, section_params in dp.summary_data.items():
            with st.expander(section, expanded=True):
                st.markdown(self.render_table(section_params), unsafe_allow_html=True)

    def samplesheet_section(self, dp):
        # Get data and show in dataframe
        samplesheet_df = dp.merged_dataframe_dict["samplesheet"]
        # samplesheet_df = samplesheet_df.reset_index(drop=True)
        st.dataframe(samplesheet_df)

    def nanopore_read_qc_section(self, dp):
        # Init
        results_dict = dp.result_dict
        dataframe_dict = dp.dataframe_dict
        # st.write("This section shows read quality reporting.")

        # Create dropdown for selecting dataset
        selected_dataset = st.selectbox("Choose a sample:", list(results_dict["toulligqc"].keys()))

        # Place charts and tables
        read_count_histogram(results_dict["toulligqc"][selected_dataset])
        read_length_scatterplot(dataframe_dict["toulligqc"][selected_dataset])

    def contaminant_removal_section(self, dp):
        # Prepare data
        host_df = dp.merged_dataframe_dict["samtools_host"]
        host_df = host_df.reset_index(drop=True)

        if "samtools_contam" in dp.merged_dataframe_dict:
            contam_df = dp.merged_dataframe_dict["samtools_contam"]
            contam_columns = contam_df.columns[1:].tolist()
            contam_columns.remove("Unmapped")

            combined_df = contam_df.reset_index(drop=True)
            combined_df.insert(2, "Host", host_df["Mapped"])
            combined_df.insert(1, "Total", contam_df.iloc[:, 1:].sum(axis=1))
            combined_columns = combined_df.columns[2:].tolist()
            combined_columns.remove("Unmapped")

        # Place charts and tables
        if "samtools_contam" in dp.merged_dataframe_dict:
            mqc_samtools_contig_bar_plot(combined_df, "Summary", combined_columns)

        mqc_samtools_bar_plot(dp.merged_dataframe_dict["samtools_host"], "Host Alignment")

        if "samtools_contam" in dp.merged_dataframe_dict:
            mqc_samtools_contig_bar_plot(contam_df, "Contaminent Alignment", contam_columns)

    def alignment_section(self, dp):
        # Place charts and tables
        mqc_samtools_bar_plot(dp.merged_dataframe_dict["samtools_align"], "Target Alignment")

    def coverage_section(self, dp):
        # Get data
        coverage_data = dp.dataframe_dict["coverage_per_base"]

        # Create dropdown for selecting dataset
        selected_dataset = st.selectbox("Choose a sample:", list(coverage_data.keys()))

        # Place chart for each contig
        for contig, df in coverage_data[selected_dataset].items():
            coverage_plot(df, contig)

    def truncation_section(self, dp):
        # Create dropdown for selecting dataset
        selected_dataset = st.selectbox("Choose a sample:", list(dp.result_dict["truncation"].keys()))

        truncation_scatterplot(dp.result_dict["truncation"][selected_dataset], "Truncation Histogram", zoom_y=False)
        truncation_scatterplot(dp.result_dict["truncation"][selected_dataset], "Truncation Histogram (Scaled)", zoom_y=True)
        truncation_barplot(dp.result_dict["truncation_type_simple"][selected_dataset], "Truncation Type Simple")
        truncation_barplot(dp.result_dict["truncation_type"][selected_dataset], "Truncation Type")

    def consensus_section(self, dp):
        # Create dropdown for selecting dataset
        selected_dataset = st.selectbox("Choose a sample:", list(dp.result_dict["consensus"].keys()))

        # Get data
        fasta_seq = "".join(dp.result_dict["consensus"][selected_dataset])
        count_table = dp.dataframe_dict["count_table"][selected_dataset]

        # Show the count table
        st.dataframe(count_table)

        # Show the fasta sequence
        st.download_button(
            label="Download Consensus FASTA",
            data=fasta_seq,
            file_name=f"{selected_dataset}.consensus.fasta",
            mime="text/plain",
        )
        st.code(f"{fasta_seq}")

    def genome_viewer_section(self, dp):
        # Create dropdown for selecting dataset
        selected_dataset = st.selectbox("Choose a sample:", list(dp.result_dict["variants_gz"].keys()))

        # Select the reference sequence and index
        if len(dp.result_dict["reference"]) > 1:
            selected_ref = selected_dataset
            selected_anno = selected_dataset
        else:
            selected_ref = list(dp.result_dict["reference"].keys())[0]
            selected_anno = list(dp.result_dict["annotation"].keys())[0]

        # Write data to static folder
        ref_path = self.tmp_dir + "_" + selected_ref + ".fasta"
        with open(ref_path, "wb") as f:
            for line in dp.result_dict["reference"][selected_ref]:
                f.write(line.encode("utf-8"))
        log.info(f"Reference written to {ref_path}")
        index_path = self.tmp_dir + "_" + selected_ref + ".fai"
        with open(index_path, "wb") as f:
            for line in dp.result_dict["reference_index"][selected_ref]:
                f.write(line.encode("utf-8"))
        log.info(f"Index written to {index_path}")
        for tool_name in dp.result_dict["variants_gz"][selected_dataset].keys():
            tool_path = self.tmp_dir + "_" + selected_dataset + "_" + tool_name + ".vcf.gz"
            with open(tool_path, "wb") as f:
                f.write(dp.result_dict["variants_gz"][selected_dataset][tool_name])
        for tool_name in dp.result_dict["variants_tbi"][selected_dataset].keys():
            tool_path = self.tmp_dir + "_" + selected_dataset + "_" + tool_name + ".vcf.gz.tbi"
            with open(tool_path, "wb") as f:
                f.write(dp.result_dict["variants_tbi"][selected_dataset][tool_name])
        ann_path = self.tmp_dir + "_" + selected_anno + ".gff"
        with open(ann_path, "wb") as f:
            for line in dp.result_dict["annotation"][selected_anno]:
                f.write(line.encode("utf-8"))
        log.info(f"Annotation written to {ann_path}")

        # Construct Uris
        base_uri = "ORIGIN_PLACEHOLDER/app/static/tmp/" + self.tmp_dir.split("/")[-1] + "_" + selected_dataset
        fasta_uri = base_uri + ".fasta"
        anno_uri = base_uri + ".gff"

        if len(dp.result_dict["reference"]) == 1:
            fasta_uri = "ORIGIN_PLACEHOLDER/app/static/tmp/" + self.tmp_dir.split("/")[-1] + "_" + selected_ref + ".fasta"
            anno_uri = "ORIGIN_PLACEHOLDER/app/static/tmp/" + self.tmp_dir.split("/")[-1] + "_" + selected_anno + ".gff"

        # Read the fasta file
        fasta_data = Fasta.read_fasta_file(ref_path)
        contigs = list(fasta_data.keys())

        # fasta_data = fasta_data[list(fasta_data.keys())[0]]
        contig_length = len(fasta_data)

        # Build the jbrowse config
        jbrowse_config = {
            "assembly": {
                "name": f"{contigs[0]}",
                "sequence": {
                    "type": "ReferenceSequenceTrack",
                    "trackId": f"{contigs[0]}-ReferenceSequenceTrack",
                    "adapter": {"type": "IndexedFastaAdapter", "uri": f"{fasta_uri}"},
                },
            },
            "defaultSession": {
                "name": "Default session",
                "margin": 0,
                "view": {
                    "id": "linearGenomeView",
                    "type": "LinearGenomeView",
                    "init": {
                        "assembly": f"{contigs[0]}",
                        "loc": f"1..{contig_length}",
                        "tracks": [
                            f"{contigs[0]}-ReferenceSequenceTrack",
                            "features",
                        ],
                    },
                },
            },
            "tracks": [
                {
                    "type": "FeatureTrack",
                    "trackId": "features",
                    "name": "Features",
                    "assemblyNames": [f"{contigs[0]}"],
                    "adapter": {
                        "type": "Gff3Adapter",
                        "uri": f"{anno_uri}",
                    },
                }
            ],
        }

        # Add variant tracks
        for tool_name in dp.result_dict["variants_gz"][selected_dataset].keys():
            jbrowse_config["tracks"].append(
                {
                    "type": "VariantTrack",
                    "trackId": f"{tool_name}_vcf_track",
                    "name": f"{tool_name} Variants",
                    "assemblyNames": [f"{contigs[0]}"],
                    "adapter": {
                        "type": "VcfTabixAdapter",
                        "uri": f"{base_uri}_{tool_name}.vcf.gz",
                    },
                }
            )

            # Also append to default session view
            jbrowse_config["defaultSession"]["view"]["init"]["tracks"].append(f"{tool_name}_vcf_track")

        # Dump into string
        config_str = json.dumps(jbrowse_config, indent=4)
        log.info(config_str)
        print(config_str)

        if self.jbrowse_component is not None:
            self.jbrowse_component(key=f"jbrowse_{selected_dataset}", config=jbrowse_config, height=1200)

    def variant_viewer_section(self, dp):
        # Get data frame
        selected_dataset = st.selectbox("Choose a sample:", list(dp.result_dict["variants"].keys()))
        df = dp.dataframe_dict["variants"][selected_dataset]
        st.dataframe(self.filter_dataframe(df))

    def filter_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Adds a UI on top of a dataframe to let viewers filter columns

        Args:
            df (pd.DataFrame): Original dataframe

        Returns:
            pd.DataFrame: Filtered dataframe
        """
        modify = st.checkbox("Add filters")

        if not modify:
            return df

        df = df.copy()

        modification_container = st.container()
        with modification_container:
            to_filter_columns = st.multiselect("Filter dataframe on", df.columns)
            for column in to_filter_columns:
                _, right = st.columns((1, 20))
                # Treat columns with < 10 unique values as categorical
                if is_categorical_dtype(df[column]) or df[column].nunique() < 10:
                    user_cat_input = right.multiselect(
                        f"Values for {column}",
                        df[column].unique(),
                        default=list(df[column].unique()),
                    )
                    df = df[df[column].isin(user_cat_input)]
                elif is_numeric_dtype(df[column]):
                    _min = float(df[column].min())
                    _max = float(df[column].max())
                    step = (_max - _min) / 100
                    user_num_input = right.slider(
                        f"Values for {column}",
                        min_value=_min,
                        max_value=_max,
                        value=(_min, _max),
                        step=step,
                    )
                    df = df[df[column].between(*user_num_input)]
                elif is_datetime64_any_dtype(df[column]):
                    user_date_input = right.date_input(
                        f"Values for {column}",
                        value=(
                            df[column].min(),
                            df[column].max(),
                        ),
                    )
                    if len(user_date_input) == 2:
                        user_date_input = tuple(map(pd.to_datetime, user_date_input))
                        start_date, end_date = user_date_input
                        df = df.loc[df[column].between(start_date, end_date)]
                else:
                    user_text_input = right.text_input(
                        f"Substring or regex in {column}",
                    )
                    if user_text_input:
                        df = df[df[column].astype(str).str.contains(user_text_input)]

        return df
