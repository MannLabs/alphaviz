import os
import torch
import plotly

import pandas as pd
import numpy as np

import plotly.graph_objects as go

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
