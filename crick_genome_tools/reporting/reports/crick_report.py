"""
Base Class for generating crick streamlit reports.
"""

import logging
import pickle

import streamlit as st

from crick_genome_tools.reporting.report_data_parser import ReportDataParser


log = logging.getLogger(__name__)


class CrickReport:
    """
    Base class for generating crick streamlit reports.
    """
    def __init__(self, report_title, data_path=None, data_obj=None, tmp_dir=None):
        self.data_path = data_path
        self.report_title = report_title
        self.tmp_dir = tmp_dir

        # Custom CSS to control sidebar width and button styling
        st.markdown(
            """
            <style>
                [data-testid="stToolbar"] {
                    display: none;
                }
                [data-testid="stSidebar"] {
                    min-width: 250px;
                    max-width: 250px;
                }
            </style>
            """,
            unsafe_allow_html=True,
        )

        # Load data once and store it in session state
        if "data_parser" not in st.session_state:
            if data_path is not None:
                # Check if data_path is a pickle file
                if data_path.endswith(".pkl"):
                    with open(data_path, "rb") as f:
                        st.session_state.data_parser = pickle.load(f)
                else:
                    # Load fresh data
                    st.session_state.data_parser = ReportDataParser(data_path)
                    st.session_state.data_parser.get_data()
            else:
                st.session_state.data_parser = data_obj

    def generate_report(self, section_headers):
        """
        Generate the report.
        """
        # Output title
        st.write("# " + self.report_title)

        # Ensure session state stores the selected section
        if "selected_section" not in st.session_state:
            st.session_state.selected_section = section_headers[0]

        # Sidebar navigation init
        st.sidebar.image("The_Francis_Crick_Institute_logo.png", width=150)
        st.sidebar.title("Navigation")
        for section in section_headers:
            if st.sidebar.button(section, key=section, help=f"Go to {section}", use_container_width=True):
                st.session_state.selected_section = section

    def activate_section(self):
        """
        Activate the selected section.
        """

