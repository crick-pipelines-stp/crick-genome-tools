"""
Helper module for generating Plotly charts for the ToulligQC report.
"""

import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objs as go

from crick_genome_tools.reporting.plotly_graph_common import line_width
from crick_genome_tools.reporting.plotly_graph_common import plotly_background_color
from crick_genome_tools.reporting.plotly_graph_common import crick_colors
from crick_genome_tools.reporting.plotly_graph_common import default_graph_layout
from crick_genome_tools.reporting.plotly_graph_common import title
from crick_genome_tools.reporting.plotly_graph_common import transparent_colors
from crick_genome_tools.reporting.plotly_graph_common import xaxis
from crick_genome_tools.reporting.plotly_graph_common import yaxis
from crick_genome_tools.reporting.plotly_graph_common import dataFrame_to_html
from crick_genome_tools.reporting.plotly_graph_common import format_int
from crick_genome_tools.reporting.plotly_graph_common import format_float


def read_count_histogram(result_dict):
    """
    Plots the histogram of count of the different types of reads:
    1D read return by Guppy
    1D pass read return by Guppy (Qscore >= 7)
    1D fail read return by Guppy (Qscore < 7)
    """

    graph_name = 'Read count histogram'
    print("HERE4")

    # Histogram with barcoded read counts
    if 'basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.count' in result_dict:

        data = {
            'All reads': result_dict['basecaller.sequencing.summary.1d.extractor.read.count'],
            'Pass reads': result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
            'Fail reads': result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"],
            'Pass barcoded reads': result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.count"],
            'Fail barcoded reads': result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.barcoded.count"]
        }

        colors = [crick_colors["all"], crick_colors["pass"], crick_colors["fail"],
                  crick_colors["barcode_pass"], crick_colors["barcode_fail"]]

        trace = go.Bar(x=[*data], y=list(data.values()),
                       hovertemplate='<b>%{x}</b><br>%{y:,}<extra></extra>',
                       marker_color=transparent_colors(colors, plotly_background_color, .5),
                       marker_line_color=colors,
                       marker_line_width=line_width)

        # Array of data for HTML table with barcode reads
        array = np.array(
            # count
            [[result_dict["basecaller.sequencing.summary.1d.extractor.read.count"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.count"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.barcoded.count"]],
             # frequencies
             [result_dict["basecaller.sequencing.summary.1d.extractor.read.count.frequency"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.frequency"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.frequency"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.frequency"],
              result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.barcoded.frequency"]]])

        dataframe = pd.DataFrame(array, index=['count', 'percent'],
                                 columns=["All reads", "Pass reads", "Fail reads", "Pass barcoded reads",
                                          "Fail barcoded reads"])

    # Histogram without barcodes
    else:

        data = {
            'All reads': result_dict['basecaller.sequencing.summary.1d.extractor.read.count'],
            'Pass reads': result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
            'Fail reads': result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"]
        }

        colors = [crick_colors['all'], crick_colors['pass'], crick_colors['fail']]

        trace = go.Bar(x=[*data], y=list(data.values()),
                       hovertemplate='<b>%{x}</b><br>%{y:,}<extra></extra>',
                       marker_color=transparent_colors(colors, plotly_background_color, .5),
                       marker_line_color=colors,
                       marker_line_width=line_width)

        # Array of data for HTML table without barcode reads
        array = np.array([[result_dict["basecaller.sequencing.summary.1d.extractor.read.count"],
                           result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
                           result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"]],
                          # frequencies
                          [result_dict["basecaller.sequencing.summary.1d.extractor.read.count.frequency"],
                           result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.frequency"],
                           result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.frequency"]]])

        # Create dataframe with array data
        dataframe = pd.DataFrame(array, index=['count', 'percent'],
                                 columns=["All reads", "Pass reads", "Fail reads"])

    # layout = go.Layout(
    #     **title(graph_name),
    #     **default_graph_layout,
    #     hovermode="x")

    layout = go.Layout(
        **title(graph_name),
        **default_graph_layout,
        hovermode="x",
        **xaxis('Read type', dict(fixedrange=True, categoryorder="total descending")),
        **yaxis('Read count'))

    fig = go.Figure(data=trace, layout=layout)
    st.plotly_chart(fig, use_container_width=True)

    # # HTML table
    # dataframe.iloc[0] = dataframe.iloc[0].astype(int).apply(lambda x: format_int(x))
    # dataframe.iloc[1:] = dataframe.iloc[1:].applymap(format_float)
    # table_html = dataFrame_to_html(dataframe)

    # return graph_name, table_html
