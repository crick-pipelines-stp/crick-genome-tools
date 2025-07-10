"""
Common functions and constants for Plotly graph generation.
"""

#  pylint: disable=invalid-name


import pkgutil
from collections import defaultdict

import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.offline as py
import streamlit as st
from scipy.interpolate import interp1d
from scipy.ndimage.filters import gaussian_filter1d


# from sklearn.utils import resample

figure_image_width = 1000
figure_image_height = 562
percent_format_str = "{:.2f}%"
line_width = 2

crick_colors = {
    "all": "#fca311",  # Yellow
    "all_1d2": "#fca311",  # Yellow
    "pass": "#51a96d",  # Green
    "fail": "#d90429",  # Red
    "barcode_pass": "#79bf90",  # Green
    "barcode_fail": "#fb1941",  # Red
    "sequence_length_over_time": "#205b47",
    "phred_score_over_time": "#7aaceb",
    "speed_over_time": "#AE3F7B",
    "nseq_over_time": "#edb773",
    "pie_chart_palette": ["#f3a683", "#f7d794", "#778beb", "#e77f67", "#cf6a87", "#786fa6", "#f8a5c2", "#63cdda", "#ea8685", "#596275"],
    "green_zone_color": "rgba(0,100,0,.1)",
}

plotly_background_color = "#e5ecf6"
legend_font_size = 16
axis_title_font_size = 14
axis_font_size = 12
on_chart_font_size = 15
title_size = 24
graph_font = "Helvetica, Arial, sans-serif"
image_dpi = 100
default_graph_layout = dict(font=dict(family=graph_font), height=figure_image_height, width=figure_image_width)
interpolation_point_count_dict = {
    "tqc_read_length_distribution": (None, 10000, 3),
    "tqc_yield_plot": (None, 10000, 3),
    "tqc_phred_score_density": (None, 1000, 3),
    "tqc_over_time_graph": (None, 1000, 3),
    "tqc_scatterplot": (10000, 4000, 3),
    "tqc_phred_violin": (10000, 4000, 3),
}


def title(graph_title):
    """
    Generates a dictionary representing the title configuration for a Plotly graph.

    Args:
        title (str): The text to be displayed as the title of the graph.

    Returns:
        dict: A dictionary containing the title configuration with text, position, and font settings.
    """
    return dict(title=dict(text=graph_title, y=0.95, x=0, xanchor="left", yanchor="top", font=dict(size=title_size, color="black")))


def transparent_component(c, b, a):
    """
    Calculate the new color component by blending the original color component with the background color component based on the alpha value.

    Args:
        c (int): The original color component (0-255).
        b (int): The background color component (0-255).
        a (float): The alpha value (0.0-1.0), where 0.0 is fully transparent and 1.0 is fully opaque.

    Returns:
        str: The new color component as a two-character hexadecimal string.
    """
    v = (1 - a) * c + a * b
    r = hex(int(v))[2:]

    if len(r) == 1:
        return "0" + r
    return r


def transparent_colors(colors, background_color, a):
    """
    Apply transparency to a list of colors by blending them with a background color based on the alpha value.

    Args:
        colors (list of str): A list of colors in hexadecimal format (e.g., '#RRGGBB').
        background_color (str): The background color in hexadecimal format (e.g., '#RRGGBB').
        a (float): The alpha value (0.0-1.0), where 0.0 is fully transparent and 1.0 is fully opaque.

    Returns:
        list of str: A list of new colors with applied transparency in hexadecimal format.
    """
    result = []

    br = int(background_color[1:3], 16)
    bg = int(background_color[3:5], 16)
    bb = int(background_color[5:7], 16)

    for c in colors:
        r = int(c[1:3], 16)
        g = int(c[3:5], 16)
        b = int(c[5:7], 16)
        new_c = "#" + transparent_component(r, br, a) + transparent_component(g, bg, a) + transparent_component(b, bb, a)
        result.append(new_c)

    return result


def xaxis(axis_title, args=None):
    """
    Create a dictionary for configuring the x-axis of a Plotly graph.

    Parameters:
    axis_title (str): The title of the x-axis.
    args (dict, optional): Additional arguments to update the x-axis configuration.

    Returns:
    dict: A dictionary containing the x-axis configuration.
    """
    axis_dict = dict(
        title="<b>" + axis_title + "</b>",
        # titlefont_size=axis_title_font_size,
        tickfont_size=axis_font_size,
    )

    if args is not None:
        axis_dict.update(dict(**args))

    return dict(xaxis=axis_dict)


def yaxis(axis_title, args=None):
    """
    Create a dictionary for configuring the y-axis of a Plotly graph.

    Parameters:
    axis_title (str): The title of the y-axis.
    args (dict, optional): Additional arguments to update the y-axis dictionary. Defaults to None.

    Returns:
    dict: A dictionary with y-axis configuration for a Plotly graph.
    """
    axis_dict = dict(
        title="<b>" + axis_title + "</b>",
        # titlefont_size=axis_title_font_size,
        tickfont_size=axis_font_size,
        fixedrange=True,
    )

    if args is not None:
        axis_dict.update(dict(**args))

    return dict(yaxis=axis_dict)


def legend(legend_title="Legend", args=None):
    legend_dict = dict(
        x=1.02,
        y=0.95,
        title_text="<b>" + legend_title + "</b>",
        title=dict(font=dict(size=legend_font_size)),
        bgcolor="white",
        bordercolor="white",
        font=dict(size=legend_font_size),
    )

    if args is not None:
        legend_dict.update(dict(**args))

    return dict(legend=legend_dict)


def dataFrame_to_html(df):
    return pd.DataFrame.to_html(df, border="")


def format_int(i):
    return "{:,d}".format(i)


def format_float(f):
    try:
        f = float(f)
        i = int(f)
        dec = f - i
        if dec == 0:
            return "{:,d}".format(i)
        else:
            return "{:,}".format(f)
    except (ValueError, TypeError):
        return "0"


def format_percent(f):
    return percent_format_str.format(f)


def make_describe_dataframe(value):
    """
    Creation of a statistics table printed with the graph in report.html
    :param value: information measured (series)
    """

    desc = value.describe()
    desc.loc["count"] = desc.loc["count"].astype(int).apply(lambda x: int(format_int(x).replace(",", "")))
    desc.iloc[1:] = desc.iloc[1:].map(lambda x: float(format_float(x).replace(",", "")))
    desc.rename({"50%": "median"}, axis="index", inplace=True)

    return desc


def interpolation_points(series, graph_name):
    count = len(series)
    threshold, npoints, sigma = interpolation_point_count_dict[graph_name]

    if threshold is not None:
        if count > threshold:
            result = npoints
        else:
            result = count
    else:
        result = npoints

    return result, sigma


def smooth_data(npoints: int, sigma: int, data, min_arg=None, max_arg=None, weights=None, density=False):
    """
    Function for smmothing data with numpy histogram function
    Returns a tuple of smooth data (ndarray)
    :param data: must be array-like data
    :param npoints: number of desired points for smoothing
    :param sigma: sigma value of the gaussian filter
    """

    if min_arg is None:
        min_arg = 0 if len(data) == 0 else np.nanmin(data)

    if max_arg is None:
        max_arg = 0 if len(data) == 0 else np.nanmax(data)

    # Compute the bin
    bins = np.linspace(min_arg, max_arg, num=npoints)

    # Compute the histogram
    y, bin_edges = np.histogram(a=data, bins=bins, weights=weights, density=density)

    # Cumulative Y
    cum_y = np.cumsum(y)

    # Center histogram
    x = bin_edges[:-1] + np.diff(bin_edges) / 2

    if min_arg == 0:
        x = np.insert(x, 0, 0)
        y = np.insert(y, 0, 0)

    if density:
        y = gaussian_filter1d(y * len(data), sigma=sigma)
        cum_y = gaussian_filter1d(cum_y * len(data), sigma=sigma)
    else:
        y = gaussian_filter1d(y, sigma=sigma)
        cum_y = gaussian_filter1d(cum_y, sigma=sigma)

    return x, y, cum_y


def read_length_distribution(graph_name, all_reads, pass_reads, fail_reads, all_color, pass_color, fail_color, xaxis_title):
    npoints, sigma = interpolation_points(all_reads, "tqc_read_length_distribution")
    min_all_reads = min(all_reads)
    max_all_reads = max(all_reads)

    count_x1, count_y1, cum_count_y1 = smooth_data(npoints=npoints, sigma=sigma, data=all_reads, min_arg=min_all_reads, max_arg=max_all_reads)
    count_x2, count_y2, cum_count_y2 = smooth_data(npoints=npoints, sigma=sigma, data=pass_reads, min_arg=min_all_reads, max_arg=max_all_reads)
    count_x3, count_y3, cum_count_y3 = smooth_data(npoints=npoints, sigma=sigma, data=fail_reads, min_arg=min_all_reads, max_arg=max_all_reads)

    sum_x1, sum_y1, cum_sum_y1 = smooth_data(
        npoints=npoints, sigma=sigma, data=all_reads, weights=all_reads, min_arg=min_all_reads, max_arg=max_all_reads
    )
    sum_x2, sum_y2, cum_sum_y2 = smooth_data(
        npoints=npoints, sigma=sigma, data=pass_reads, weights=pass_reads, min_arg=min_all_reads, max_arg=max_all_reads
    )
    sum_x3, sum_y3, cum_sum_y3 = smooth_data(
        npoints=npoints, sigma=sigma, data=fail_reads, weights=fail_reads, min_arg=min_all_reads, max_arg=max_all_reads
    )

    # Find 50 percentile for zoomed range on x axis
    max_x_range = np.percentile(all_reads, 99)

    coef = max_all_reads / npoints

    max_y = max(max(count_y1), max(count_y2), max(count_y3)) / coef
    max_sum_y = max(max(sum_y1), max(sum_y2), max(sum_y3)) / coef

    fig = go.Figure()

    # Read graphs
    fig.add_trace(go.Scatter(x=count_x1, y=count_y1 / coef, name="All reads", fill="tozeroy", marker_color=all_color, visible=True))
    fig.add_trace(go.Scatter(x=count_x2, y=count_y2 / coef, name="Pass reads", fill="tozeroy", marker_color=pass_color, visible=True))
    fig.add_trace(go.Scatter(x=count_x3, y=count_y3 / coef, name="Fail reads", fill="tozeroy", marker_color=fail_color, visible=True))

    # Threshold
    for p in [25, 50, 75]:
        x0 = np.percentile(all_reads, p)
        if p == 50:
            t = "median<br>all reads"
        else:
            t = str(p) + "%<br>all reads"
        fig.add_trace(
            go.Scatter(
                mode="lines+text",
                name="All reads",
                x=[x0, x0],
                y=[0, max_y],
                line=dict(color="gray", width=1, dash="dot"),
                text=["", t],
                textposition="top center",
                hoverinfo="skip",
                showlegend=False,
                visible=True,
            )
        )

    # Base plots
    # Read graphs
    fig.add_trace(go.Scatter(x=sum_x1, y=sum_y1 / coef, name="All reads", fill="tozeroy", marker_color=all_color, visible=False))
    fig.add_trace(go.Scatter(x=sum_x2, y=sum_y2 / coef, name="Pass reads", fill="tozeroy", marker_color=pass_color, visible=False))
    fig.add_trace(go.Scatter(x=sum_x3, y=sum_y3 / coef, name="Fail reads", fill="tozeroy", marker_color=fail_color, visible=False))

    # Threshold
    for p in [25, 50, 75]:
        x0 = np.percentile(all_reads, p)
        if p == 50:
            t = "median<br>all reads"
        else:
            t = str(p) + "%<br>all reads"
        fig.add_trace(
            go.Scatter(
                mode="lines+text",
                name="All reads",
                x=[x0, x0],
                y=[0, max_y],
                line=dict(color="gray", width=1, dash="dot"),
                text=["", t],
                textposition="top center",
                hoverinfo="skip",
                showlegend=False,
                visible=False,
            )
        )

    fig.update_layout(
        **title(graph_name),
        **default_graph_layout,
        **legend(args=dict(y=0.75)),
        hovermode="x",
        **xaxis(xaxis_title, dict(range=[min_all_reads, max_x_range], type="linear")),
        **yaxis("Read count", dict(range=[0, max_y * 1.10])),
    )

    # Add buttons
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                direction="down",
                buttons=list(
                    [
                        dict(
                            args=[
                                {"visible": [True, True, True, True, True, True, False, False, False]},
                                {
                                    "xaxis": {"type": "linear", "range": [min_all_reads, max_x_range]},
                                    "yaxis": {"title": "<b>Read count</b>", "range": [0, max_y * 1.10]},
                                },
                            ],
                            label="Reads linear",
                            method="update",
                        ),
                        dict(
                            args=[
                                {"visible": [True, True, True, True, True, True, False, False, False]},
                                {"xaxis": {"type": "log"}, "yaxis": {"title": "<b>Read count</b>", "range": [0, max_y * 1.10]}},
                            ],
                            label="Reads log",
                            method="update",
                        ),
                        dict(
                            args=[
                                {"visible": [False, False, False, False, False, False, True, True, True]},
                                {
                                    "xaxis": {"type": "linear", "range": [min_all_reads, max_x_range]},
                                    "yaxis": {"title": "<b>Base count</b>", "range": [0, max_sum_y * 1.10]},
                                },
                            ],
                            label="Bases linear",
                            method="update",
                        ),
                        dict(
                            args=[
                                {"visible": [False, False, False, False, False, False, True, True, True]},
                                {"xaxis": {"type": "log"}, "yaxis": {"title": "<b>Base count</b>", "range": [0, max_sum_y * 1.10]}},
                            ],
                            label="Bases log",
                            method="update",
                        ),
                    ]
                ),
                pad={"r": 20, "t": 20, "l": 20, "b": 20},
                showactive=True,
                x=1.0,
                xanchor="left",
                y=1.25,
                yanchor="top",
            ),
        ]
    )
    st.plotly_chart(fig, use_container_width=True)

    table_df = pd.concat([pd.Series(all_reads), pass_reads, fail_reads], axis=1, keys=["All reads", "Pass reads", "Fail reads"])
    table_df = make_describe_dataframe(table_df)
    st.dataframe(table_df, use_container_width=True)
