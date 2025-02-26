"""
Base Class for generating crick streamlit reports.
"""

import logging
from crick_genome_tools.reporting.report_data_parser import ReportDataParser

import streamlit as st

log = logging.getLogger(__name__)

class CrickReport:
    def __init__(self, data_path):
        self.data_path = data_path
        self.data_parser = ReportDataParser(data_path)
        self.data_parser.get_data()

    def generate_report(self, report_title):
        """
        Generate the report.
        """

        # Set title and favicon
        st.set_page_config(
            page_title=report_title,
            page_icon="favicon.ico",
        )
        st.write("# " + report_title)
