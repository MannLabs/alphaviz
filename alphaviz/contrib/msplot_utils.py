import os
import plotly

import pandas as pd
import numpy as np

import plotly.graph_objects as go
import plotly.express as px

def _plot_scatter(
    fig:go.Figure,
    x_values, y_values,
    color, marker_size,
    hovertext, 
    hovertemplate,
    name='',
    row=None, col=None,
):
    fig.add_trace(
        go.Scatter(
            x=x_values,
            y=y_values,
            hovertext=hovertext,
            hovertemplate=hovertemplate,
            mode='markers',
            marker=dict(
                color=color,
                size=marker_size
            ),
            name=name,
            showlegend=False
        ),
        row=row, col=col
    )


def _plot_line(
    tims_data,
    selected_indices: np.ndarray,
    label: str,
    marker_color: str,
    remove_zeros: bool = False,
    trim: bool = True,
):
    """Plot an XIC as a lineplot.

    Parameters
    ----------
    tims_data : alphatims.bruker.TimsTOF
        An alphatims.bruker.TimsTOF data object.
    selected_indices : np.ndarray
        The raw indices of tims_data that are selected for this plot.
    label : str
        The label for the line plot.
    remove_zeros : bool
        If True, zeros are removed. Default: False.
    trim : bool
        If True, zeros on the left and right are trimmed. Default: True.

    Returns
    -------
    go.Figure
        the XIC line plot.
    """
    axis_dict = {
        "rt": "RT (min)",
        "intensity": "Intensity",
    }
    x_axis_label = axis_dict["rt"]
    # y_axis_label = axis_dict["intensity"]
    labels = {
        'RT (min)': "rt_values",
    }
    x_dimension = labels[x_axis_label]
    intensities = tims_data.bin_intensities(selected_indices, [x_dimension])
    x_ticks = tims_data.rt_values / 60

    non_zeros = np.flatnonzero(intensities)
    if len(non_zeros) == 0:
        x_ticks = np.empty(0, dtype=x_ticks.dtype)
        intensities = np.empty(0, dtype=intensities.dtype)
    else:
        if remove_zeros:
            x_ticks = x_ticks[non_zeros]
            intensities = intensities[non_zeros]
        elif trim:
            start = max(0, non_zeros[0] - 1)
            end = non_zeros[-1] + 2
            x_ticks = x_ticks[start: end]
            intensities = intensities[start: end]

    trace = go.Scatter(
        x=x_ticks,
        y=intensities,
        mode='lines',
        text=[f'{x_axis_label}'.format(i + 1) for i in range(len(x_ticks))],
        hovertemplate='<b>%{text}:</b> %{x};<br><b>Intensity:</b> %{y}',
        name=label,
        marker=marker_color,
        legendgroup=label.split(' ')[0],
    )
    return trace

