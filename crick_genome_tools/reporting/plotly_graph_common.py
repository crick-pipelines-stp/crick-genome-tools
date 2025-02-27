"""
Common functions and constants for Plotly graph generation.
"""

#  pylint: disable=invalid-name


from collections import defaultdict

import pkgutil
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.offline as py

figure_image_width = 1000
figure_image_height = 562
percent_format_str = '{:.2f}%'
line_width = 2

crick_colors = {'all': '#fca311',  # Yellow
                    'all_1d2': '#fca311',  # Yellow
                    'pass': '#51a96d',  # Green
                    'fail': '#d90429',  # Red
                    'barcode_pass': '#79bf90',  # Green
                    'barcode_fail': '#fb1941',  # Red
                    'sequence_length_over_time': '#205b47',
                    'phred_score_over_time': '#7aaceb',
                    'speed_over_time': '#AE3F7B',
                    'nseq_over_time': '#edb773',
                    'pie_chart_palette': ["#f3a683", "#f7d794", "#778beb", "#e77f67", "#cf6a87", "#786fa6", "#f8a5c2",
                                          "#63cdda", "#ea8685", "#596275"],
                    'green_zone_color': 'rgba(0,100,0,.1)'
                    }

plotly_background_color = '#e5ecf6'
legend_font_size = 16
axis_title_font_size = 14
axis_font_size = 12
on_chart_font_size = 15
title_size = 24
graph_font = 'Helvetica, Arial, sans-serif'
image_dpi = 100
default_graph_layout = dict(
    font=dict(family=graph_font),
    height=figure_image_height,
    width=figure_image_width
)

def title(graph_title):
    """
    Generates a dictionary representing the title configuration for a Plotly graph.

    Args:
        title (str): The text to be displayed as the title of the graph.

    Returns:
        dict: A dictionary containing the title configuration with text, position, and font settings.
    """
    return dict(title=dict(
        text=graph_title,
        y=0.95,
        x=0,
        xanchor='left',
        yanchor='top',
        font=dict(size=title_size, color="black")))


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
        return '0' + r
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
        new_c = '#' + \
                transparent_component(r, br, a) + \
                transparent_component(g, bg, a) + \
                transparent_component(b, bb, a)
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
        title='<b>' + axis_title + '</b>',
        titlefont_size=axis_title_font_size,
        tickfont_size=axis_font_size)

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
        title='<b>' + axis_title + '</b>',
        titlefont_size=axis_title_font_size,
        tickfont_size=axis_font_size,
        fixedrange=True)

    if args is not None:
        axis_dict.update(dict(**args))

    return dict(yaxis=axis_dict)


def dataFrame_to_html(df):
    return pd.DataFrame.to_html(df, border="")


def format_int(i):
    return '{:,d}'.format(i)


def format_float(f):
    try:
        s = str(f)
        i = int(s.split('.')[0])
        f = float('0.' + s.split('.')[1])
    except:
        return 0

    return '{:,d}'.format(i) + '{:.2f}'.format(f)[1:]


def format_percent(f):
    return percent_format_str.format(f)