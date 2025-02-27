"""
Class for generating vector core AAV report.
"""

import logging
from crick_genome_tools.reporting.reports.crick_report import CrickReport

import streamlit as st

log = logging.getLogger(__name__)

class VectorCoreAavReport(CrickReport):

    def __init__(self, data_path):
        super().__init__(data_path, "Vectorcore AAV Report")

    def generate_report(self, section_headers):
        section_headers = [
            "Introduction",
            "Data Overview",
            "Analysis",
            "Results",
            "Conclusion"
        ]
        super().generate_report(section_headers)

        # Activate current section
        self.activate_section()

    def activate_section(self):
        # Display the selected section content
        st.title(st.session_state.selected_section)

        if st.session_state.selected_section == "Introduction":
            st.write("Welcome to the report. This section provides an overview of the study.")

        elif st.session_state.selected_section == "Data Overview":
            st.write("This section covers data sources, preprocessing, and key statistics.")

        elif st.session_state.selected_section == "Analysis":
            st.write("Here, we dive into the data analysis with visualizations and insights.")

        elif st.session_state.selected_section == "Results":
            st.write("Summary of the findings and key takeaways.")

        elif st.session_state.selected_section == "Conclusion":
            st.write("Final thoughts and recommendations based on the analysis.")
