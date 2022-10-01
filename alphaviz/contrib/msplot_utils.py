import os
import plotly

import pandas as pd
import numpy as np

import plotly.graph_objects as go
import plotly.express as px

from alphatims.bruker import TimsTOF

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

def _plot_line_indices(
    tims_data:TimsTOF,
    selected_indices: np.ndarray,
    label: str,
    marker_color: str,
    remove_zeros: bool = False,
    trim: bool = True,
    view_dim: str='rt' # or 'im'
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
    labels = {
        'rt': "rt_values",
        'im': "mobility_values"
    }

    x_dimension = labels[view_dim]

    intensities = tims_data.bin_intensities(selected_indices, [x_dimension])
    if view_dim == 'rt': 
        x_ticks = tims_data.rt_values / 60
    else: 
        x_ticks = tims_data.mobility_values

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
        text=[
            f'RT: {_x*60:.3f}s' for _x in x_ticks] if view_dim == 'rt' 
            else [f'IM: {_x:.3f}' for _x in x_ticks],
        hovertemplate='%{text} <br><b>Intensity:</b> %{y}',
        name=label,
        marker=marker_color,
        legendgroup=label.split(' ')[0],
    )
    return trace

def _plot_line(
    tims_sliced_df:pd.DataFrame,
    frame_df:pd.DataFrame,
    label: str,
    marker_color: str,
    view_dim: str='rt' # or 'im'
):
    """Plot an XIC as a lineplot.

    Parameters
    ----------
    tims_sliced_df : pd.DataFrame
        TimsTOF[...] df
    label : str
        The label for the line plot.
    trim : bool
        If True, zeros on the left and right are trimmed. Default: True.

    Returns
    -------
    go.Figure
        the XIC line plot.
    """

    if view_dim == 'rt':
        tims_sliced_df = tims_sliced_df.groupby(
            'frame_indices', as_index=False
        ).agg(
            {
                'rt_values':'mean',
                'intensity_values':'sum',
            }
        )
        tims_sliced_df['rt_values'] /= 60
        tims_sliced_df.sort_values('rt_values', inplace=True)
        tims_sliced_df = frame_df.merge(tims_sliced_df, on=['frame_indices','rt_values'], how='left')
        tims_sliced_df.loc[tims_sliced_df.intensity_values.isna(),'intensity_values'] = 0
        x_ticks = tims_sliced_df.rt_values.values
    else: 
        tims_sliced_df = tims_sliced_df.groupby(
            'scan_indices', as_index=False
        ).agg(
            {
                'mobility_values':'mean',
                'intensity_values':'sum',
            }
        )
        tims_sliced_df.sort_values('mobility_values', inplace=True)
        x_ticks = tims_sliced_df.mobility_values.values

    trace = go.Scatter(
        x=x_ticks,
        y=tims_sliced_df.intensity_values.values,
        mode='lines',
        text=[
            f'RT: {_x*60:.3f}s' for _x in x_ticks] if view_dim == 'rt' 
            else [f'IM: {_x:.3f}' for _x in x_ticks],
        hovertemplate='%{text}<br>Intensity: %{y}',
        name=label,
        marker=marker_color,
        legendgroup=label.split(' ')[0],
    )
    return trace