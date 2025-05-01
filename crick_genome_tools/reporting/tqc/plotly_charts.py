"""
Helper module for generating Plotly charts.
"""

import numpy as np
import pandas as pd
import plotly.graph_objs as go
import streamlit as st

from crick_genome_tools.reporting.plotly_graph_common import (
    crick_colors,
    default_graph_layout,
    format_float,
    format_int,
    line_width,
    plotly_background_color,
    read_length_distribution,
    title,
    transparent_colors,
    xaxis,
    yaxis,
    smooth_data,
    legend
)


def read_count_histogram(result_dict):
    """
    Plots the histogram of count of the different types of reads:
    1D read return by Guppy
    1D pass read return by Guppy (Qscore >= 7)
    1D fail read return by Guppy (Qscore < 7)
    """

    graph_name = "Read count histogram"

    # Histogram with barcoded read counts
    if "basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.count" in result_dict:

        data = {
            "All reads": result_dict["basecaller.sequencing.summary.1d.extractor.read.count"],
            "Pass reads": result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
            "Fail reads": result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"],
            "Pass barcoded reads": result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.count"],
            "Fail barcoded reads": result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.barcoded.count"],
        }

        colors = [crick_colors["all"], crick_colors["pass"], crick_colors["fail"], crick_colors["barcode_pass"], crick_colors["barcode_fail"]]

        trace = go.Bar(
            x=[*data],
            y=list(data.values()),
            hovertemplate="<b>%{x}</b><br>%{y:,}<extra></extra>",
            marker_color=transparent_colors(colors, plotly_background_color, 0.5),
            marker_line_color=colors,
            marker_line_width=line_width,
        )

        # Array of data for HTML table with barcode reads
        array = np.array(
            # count
            [
                [
                    result_dict["basecaller.sequencing.summary.1d.extractor.read.count"],
                    result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
                    result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"],
                    result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.count"],
                    result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.barcoded.count"],
                ],
                # frequencies
                [
                    result_dict["basecaller.sequencing.summary.1d.extractor.read.count.frequency"],
                    result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.frequency"],
                    result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.frequency"],
                    result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.barcoded.frequency"],
                    result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.barcoded.frequency"],
                ],
            ]
        )

        dataframe = pd.DataFrame(
            array, index=["count", "percent"], columns=["All reads", "Pass reads", "Fail reads", "Pass barcoded reads", "Fail barcoded reads"]
        )

    # Histogram without barcodes
    else:

        data = {
            "All reads": result_dict["basecaller.sequencing.summary.1d.extractor.read.count"],
            "Pass reads": result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
            "Fail reads": result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"],
        }

        colors = [crick_colors["all"], crick_colors["pass"], crick_colors["fail"]]

        trace = go.Bar(
            x=[*data],
            y=list(data.values()),
            hovertemplate="<b>%{x}</b><br>%{y:,}<extra></extra>",
            marker_color=transparent_colors(colors, plotly_background_color, 0.5),
            marker_line_color=colors,
            marker_line_width=line_width,
        )

        # Array of data for HTML table without barcode reads
        array = np.array(
            [
                [
                    result_dict["basecaller.sequencing.summary.1d.extractor.read.count"],
                    result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.count"],
                    result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.count"],
                ],
                # frequencies
                [
                    result_dict["basecaller.sequencing.summary.1d.extractor.read.count.frequency"],
                    result_dict["basecaller.sequencing.summary.1d.extractor.read.pass.frequency"],
                    result_dict["basecaller.sequencing.summary.1d.extractor.read.fail.frequency"],
                ],
            ]
        )

        # Create dataframe with array data
        dataframe = pd.DataFrame(array, index=["count", "percent"], columns=["All reads", "Pass reads", "Fail reads"])

    layout = go.Layout(
        **title(graph_name),
        **default_graph_layout,
        hovermode="x",
        **xaxis("Read type", dict(fixedrange=True, categoryorder="total descending")),
        **yaxis("Read count"),
    )

    fig = go.Figure(data=trace, layout=layout)
    st.plotly_chart(fig, use_container_width=True)

    # HTML table
    dataframe.iloc[0] = dataframe.iloc[0].astype(int).apply(lambda x: int(format_int(x).replace(",", "")))
    dataframe.iloc[1:] = dataframe.iloc[1:].map(lambda x: float(format_float(x).replace(",", "")))
    st.dataframe(dataframe, use_container_width=True)


def read_length_scatterplot(dataframe_dict):
    """
    Generates a scatter plot for the distribution of read lengths.

    Parameters:
    dataframe_dict (dict): A dictionary containing the following keys:
        - 'all.reads.sequence.length': Data for all reads.
        - 'pass.reads.sequence.length': Data for reads that passed quality control.
        - 'fail.reads.sequence.length': Data for reads that failed quality control.

    Returns:
    plotly.graph_objs._figure.Figure: A Plotly figure object representing the scatter plot.

    """
    graph_name = "Distribution of read lengths"

    return read_length_distribution(
        graph_name=graph_name,
        all_reads=dataframe_dict["all.reads.sequence.length"],
        pass_reads=dataframe_dict["pass.reads.sequence.length"],
        fail_reads=dataframe_dict["fail.reads.sequence.length"],
        all_color=crick_colors["all"],
        pass_color=crick_colors["pass"],
        fail_color=crick_colors["fail"],
        xaxis_title="Read length (bp)",
    )


def mqc_samtools_bar_plot(data_frame, graph_name):
    # Init
    mapped_colors = [crick_colors["pass"]] * len(data_frame)
    unmapped_colors = [crick_colors["fail"]] * len(data_frame)

    trace_mapped = go.Bar(
        x=data_frame["Sample"],
        y=data_frame["Mapped"],
        hovertemplate="<b>%{x}</b><br>%{y:,}<extra></extra>",
        marker_color=transparent_colors(mapped_colors, plotly_background_color, 0.5),
        marker_line_color=mapped_colors,
        marker_line_width=line_width,
        name="Mapped",
    )

    trace_unmapped = go.Bar(
        x=data_frame["Sample"],
        y=data_frame["Unmapped"],
        hovertemplate="<b>%{x}</b><br>%{y:,}<extra></extra>",
        marker_color=transparent_colors(unmapped_colors, plotly_background_color, 0.5),
        marker_line_color=unmapped_colors,
        marker_line_width=line_width,
        name="Unmapped",
    )

    layout = go.Layout(
        **title(graph_name), **default_graph_layout, barmode="stack", hovermode="x", **xaxis("Sample", dict(fixedrange=True)), **yaxis("Read count")
    )

    fig = go.Figure(data=[trace_mapped, trace_unmapped], layout=layout)
    st.plotly_chart(fig, use_container_width=True)
    data_frame.index = data_frame["Sample"]
    data_frame = data_frame.drop(columns=["Sample"])
    st.dataframe(data_frame, use_container_width=True)


def mqc_samtools_contig_bar_plot(data_frame, graph_name, contigs):
    # Init
    mapped_colors = [crick_colors["pass"]] * len(data_frame)
    unmapped_colors = [crick_colors["fail"]] * len(data_frame)
    extra_colors = crick_colors["pie_chart_palette"]
    extra_color_idx = 0
    traces = []

    trace_mapped = go.Bar(
        x=data_frame["Sample"],
        y=data_frame[contigs[0]],
        hovertemplate="<b>%{x} - Primary</b><br>%{y:,}<extra></extra>",
        marker_color=transparent_colors(mapped_colors, plotly_background_color, 0.5),
        marker_line_color=mapped_colors,
        marker_line_width=line_width,
        name=contigs[0],
    )
    traces.append(trace_mapped)

    for contig in contigs[1:]:
        contig_colours = [extra_colors[extra_color_idx]] * len(data_frame)
        trace_extra = go.Bar(
            x=data_frame["Sample"],
            y=data_frame[contig],
            hovertemplate=f"<b>%{{x}} - {contig}</b><br>%{{y:,}}<extra></extra>",
            marker_color=transparent_colors(contig_colours, plotly_background_color, 0.5),
            marker_line_color=contig_colours,
            marker_line_width=line_width,
            name=contig,
        )
        extra_color_idx += 1
        if extra_color_idx == len(extra_colors):
            extra_color_idx = 0
        traces.append(trace_extra)

    trace_unmapped = go.Bar(
        x=data_frame["Sample"],
        y=data_frame["Unmapped"],
        hovertemplate="<b>%{x} - Unmapped</b><br>%{y:,}<extra></extra>",
        marker_color=transparent_colors(unmapped_colors, plotly_background_color, 0.5),
        marker_line_color=unmapped_colors,
        marker_line_width=line_width,
        name="Unmapped",
    )
    traces.append(trace_unmapped)

    layout = go.Layout(
        **title(graph_name), **default_graph_layout, barmode="stack", hovermode="x", **xaxis("Sample", dict(fixedrange=True)), **yaxis("Read count")
    )

    fig = go.Figure(data=traces, layout=layout)
    st.plotly_chart(fig, use_container_width=True)
    data_frame.index = data_frame["Sample"]
    data_frame = data_frame.drop(columns=["Sample"])
    st.dataframe(data_frame, use_container_width=True)


def coverage_plot(data_frame, graph_name):
    trace = go.Scatter(
        x=data_frame["Position"], y=data_frame["Depth"], mode="lines", hovertemplate="<b>Position</b>: %{x}<br><b>Depth</b>: %{y:,}<extra></extra>"
    )

    layout = go.Layout(**title(graph_name), **default_graph_layout, hovermode="x", **xaxis("Position", dict(fixedrange=True)), **yaxis("Depth"))

    fig = go.Figure(data=trace, layout=layout)
    st.plotly_chart(fig, use_container_width=True)


def truncation_scatterplot(data_frame):
    # Init
    graph_name = "Truncation histogram"
    npoints = 10000
    sigma = 3
    start_pos = data_frame["Read Start"]
    end_pos = data_frame["Read End"]
    min_all_pos = 0
    max_all_pos = max(start_pos.max(), end_pos.max())

    count_x1, count_y1, cum_count_y1 = smooth_data(npoints=npoints, sigma=sigma, data=start_pos, min_arg=min_all_pos, max_arg=max_all_pos)
    count_x2, count_y2, cum_count_y2 = smooth_data(npoints=npoints, sigma=sigma, data=end_pos, min_arg=min_all_pos, max_arg=max_all_pos)
    max_y = max(max(count_y1), max(count_y2))

    # Calc upper y
    combined_y = np.concatenate([count_y1, count_y2])
    y_upper = np.percentile(combined_y, 99.5)

    fig = go.Figure()

    # Read graphs
    fig.add_trace(go.Scatter(x=count_x1, y=count_y1, name="Read Start", fill="tozeroy", marker_color=crick_colors["pass"], visible=True))
    fig.add_trace(go.Scatter(x=count_x2, y=count_y2, name="Read End", fill="tozeroy", marker_color=crick_colors["fail"], visible=True))

    fig.update_layout(
        **title(graph_name),
        **default_graph_layout,
        **legend(args=dict(y=0.75)),
        hovermode="x",
        **xaxis("Reference Position", dict(range=[min_all_pos, max_all_pos], type="linear")),
        **yaxis("Read count", dict(range=[0, y_upper])),
    )

    st.plotly_chart(fig, use_container_width=True)


def truncation_barplot(data_frame):
    # Init
    graph_name = "Truncation Type"
    counts = data_frame["aln_type"].value_counts().sort_values(ascending=False)
    data_frame = counts.reset_index()
    data_frame.columns = ["Truncation Type", "Count"]
    total = counts.sum()
    data_frame["Percentage (%)"] = (data_frame["Count"] / total * 100).round(2)
    colors = [crick_colors["pie_chart_palette"][2]] * len(data_frame)

    trace = go.Bar(
        x=data_frame["Truncation Type"],
        y=data_frame["Count"],
        hovertemplate="<b>%{x}</b><br>%{y:,} reads<extra></extra>",
        marker_color=transparent_colors(colors, plotly_background_color, 0.5),
        marker_line_color=colors,
        marker_line_width=line_width,
        name="Alignment Type"
    )

    layout = go.Layout(
        **title(graph_name),
        **default_graph_layout,
        barmode="stack",
        hovermode="x",
        **xaxis("Truncation type", dict(tickangle=45, fixedrange=True)),
        **yaxis("Read count")
    )

    fig = go.Figure(data=[trace], layout=layout)
    st.plotly_chart(fig, use_container_width=True)
    st.dataframe(data_frame, use_container_width=True, hide_index=True)
