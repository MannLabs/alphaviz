#!python
"""
This module provides the plotting functions used by AlphaViz.
"""
import re

import numpy as np
import pandas as pd

import plotly.graph_objects as go
import plotly.subplots
import plotly.express as px

import holoviews as hv
from holoviews.operation.datashader import dynspread, rasterize, shade
from bokeh.io import export_svgs

import alphaviz.preprocessing
import alphaviz.utils


def plot_sequence_coverage(
    sequence: str,
    gene_name: str,
    peptides_list: list,
    colorscale_qualitative: str,
    colorscale_sequential: str,
    regex: str,
    prot_id: str = ""
) -> go.Figure:
    """Create a protein sequence coverage plot.

    Parameters
    ----------
    sequence : str
        Protein sequence.
    gene_name : str
        Gene name.
    peptides_list : list
        List of all identified peptides for the specified protein.
    colorscale_qualitative : str
        A name of a built-in qualitative Plotly color scale.
    colorscale_sequential : str
        A name of a built-in sequential Plotly color scale.
    regex : str
        A regular expression to be applied for the peptide sequence.
    Returns
    -------
    plotly.graph_objects.Figure object
        A protein sequence coverage plot showing all peptides on the protein sequence.

    """
    fig = go.Figure()

    try:
        aa_position = list(range(1, len(sequence) + 1))
        aa_list = list(sequence)
    except TypeError:
        print(f"The AA sequence is not found in the fasta file for the selected protein_id {prot_id}.")
        return None

    fig.add_trace(
        go.Bar(
            x=aa_position,
            y=np.ones(len(aa_position))*0.5,
            name='Protein aa sequence',
            hovertext=aa_list,
            hovertemplate='<b>amino acid:</b> %{hovertext};<br><b>position:</b> %{x}.',
            marker_color='grey'
        )
    )
    selected_peptide_cov = np.zeros(len(sequence), dtype=np.bool)
    if len(peptides_list) <= len(getattr(px.colors.qualitative, colorscale_qualitative)):
        for ind, peptide_mod in enumerate(peptides_list):
            peptide_mod = peptide_mod.replace('_', '')
            peptide = re.sub(regex, "", peptide_mod)
            start = sequence.find(peptide)
            peptide_cov = range(start + 1, start + len(peptide) + 1)
            selected_peptide_cov[start + 1: start + len(peptide) + 1] = True
            if start != -1:
                fig.add_trace(
                    go.Bar(
                        x=list(peptide_cov),
                        y=np.ones(len(peptide))*2,
                        name=f'Peptide: {peptide_mod}',
                        hovertemplate='<br><b>position:</b> %{x}.',
                        opacity=0.5,
                        marker=dict(color=getattr(px.colors.qualitative, colorscale_qualitative)[ind])
                    )
                )
            else:
                print(f'The peptide {peptide} is not found in the protein sequence of the protein_id {prot_id}.')
                return None
    else:
        colorscale_sequential_colors = px.colors.sample_colorscale(colorscale_sequential, samplepoints=len(peptides_list))
        for ind, peptide_mod in enumerate(peptides_list):
            peptide_mod = peptide_mod.replace('_', '')
            peptide = re.sub(regex, "", peptide_mod)
            start = sequence.find(peptide)
            peptide_cov = range(start + 1, start + len(peptide) + 1)
            selected_peptide_cov[start + 1: start + len(peptide) + 1] = True
            if start != -1:
                fig.add_trace(
                    go.Bar(
                        x=list(peptide_cov),
                        y=np.ones(len(peptide))*2,
                        name=f'Peptide: {peptide_mod}',
                        hovertemplate='<br><b>position:</b> %{x}.',
                        opacity=0.5,
                        marker=dict(color=colorscale_sequential_colors[ind])
                    )
                )
            else:
                print(f'The peptide {peptide} is not found in the protein sequence of the protein_id {prot_id}.')
                return None
    aa_coverage = round(np.sum(selected_peptide_cov) / len(selected_peptide_cov) * 100, 2)
    fig.update_layout(
        title=dict(
            text=f"Protein coverage diagram (protein ID {prot_id})" + '\n' + f"(AA coverage {aa_coverage}%)",
            font=dict(
                size=16,
            ),
            y=0.99,
            x=0.5,
            xanchor='center',
            yanchor='top'
        ),
        yaxis=dict(
            tickfont_size=1,
            showticklabels=False,
            nticks=0,
            visible=False
        ),
        barmode='overlay',
        bargap=0,  # gap between bars of adjacent location coordinates.
        bargroupgap=0,  # gap between bars of the same location coordinate.
        hovermode="x",
        template="plotly_white",  # "plotly", "plotly_white", "plotly_dark", "ggplot2", "seaborn", "simple_white"
        # width=1000,
        height=200
    )

    fig.update_yaxes(range=[-1, 3])
    fig.update_layout(showlegend=False)

    return fig


def plot_chrom(
    data,  # alphatims.bruker.TimsTOF object
    colorscale_qualitative: str,
) -> go.Figure:
    """Create a plot showing 4 chromatogram types: total ion chromatogram (TIC) and base peak chromatogram (BPC) for MS1 and MS2 raw data.

    Parameters
    ----------
    data : alphatims.bruker.TimsTOF object
        An Alphatims.bruker.TimsTOF object.

    Returns
    -------
    plotly.graph_objects.Figure object
        A plotly.graph_objects.Figure object showing overlapping TIC and BPC for MS1 and MS2 raw data.

    """

    chrom_ms1 = data.frames.query('MsMsType == 0')[['Time', 'SummedIntensities', 'MaxIntensity']]
    chrom_ms1['RT'] = chrom_ms1['Time'] / 60
    chrom_ms2 = data.frames.query('MsMsType != 0')[['Time', 'SummedIntensities', 'MaxIntensity']]
    chrom_ms2['RT'] = chrom_ms2['Time'] / 60

    fig = go.Figure()

    total_ion_col = ['RT', 'SummedIntensities']
    base_peak_col = ['RT', 'MaxIntensity']

    for i, chrom_type in enumerate(['Total Ion Chromatogram - MS1', 'Base Peak Chromatogram - MS1', 'Total Ion Chromatogram - MS2', 'Base Peak Chromatogram - MS2']):
        if chrom_type == 'Total Ion Chromatogram - MS1':
            data = chrom_ms1[total_ion_col]
        elif chrom_type == 'Total Ion Chromatogram - MS2':
            data = chrom_ms2[total_ion_col]
        elif chrom_type == 'Base Peak Chromatogram - MS1':
            data = chrom_ms1[base_peak_col]
        elif chrom_type == 'Base Peak Chromatogram - MS2':
            data = chrom_ms2[base_peak_col]
        fig.add_trace(
            go.Scatter(
                x=data.iloc[:, 0],
                y=data.iloc[:, 1],
                name=chrom_type,
                hovertemplate='<b>RT:</b> %{x};<br><b>Intensity:</b> %{y}.',
                marker=dict(color=getattr(px.colors.qualitative, colorscale_qualitative)[i]),
            )
        )

    fig.update_layout(
        title=dict(
            text="Chromatograms",
            font=dict(
                size=16,
            ),
            x=0.5,
            xanchor='center',
            yanchor='top',
            y=0.92
        ),
        xaxis=dict(
            title='RT, min',
            titlefont_size=14,
            tickmode='auto',
            tickfont_size=14,
            autorange=False,
        ),
        yaxis=dict(
            title='Intensity',
            # type="log"
        ),
        legend=dict(
            orientation="h",
            x=1,
            xanchor="right",
            yanchor="bottom",
            y=1.01
        ),
        legend_title_text='Types:',
        hovermode="x",
        template="plotly_white",
        # width=1000,
        height=450
    )

    fig.update_xaxes(range=[0, max(chrom_ms1.RT.max(), chrom_ms2.RT.max())])
    return fig


def plot_heatmap(
    df: pd.DataFrame,
    x_axis_label: str = "m/z, Th",
    y_axis_label: str = "Inversed IM, V·s·cm\u207B\u00B2",
    z_axis_label: str = "Intensity",
    mz: float = 0.0,
    im: float = 0.0,
    title: str = "",
    width: int = 450,
    height: int = 450,
    background_color: str = "black",
    precursor_color: str = 'green',
    precursor_size: int = 15,
    colormap: str = 'fire',
    **kwargs
) -> hv.Scatter:
    """Create a heatmap for the MS1/MS2 frame that overlaps with the precursor mark at the location where the precursor has been selected for analysis.

    Parameters
    ----------
    df : pd.DataFrame
        A slice of the alphatims.bruker.TimsTOF object with MS1 or MS2/PASEF frame data.
    mz: float
        M/z value of the precursor identified in the frame.
    im: float
        IM value of the precursor identified in the frame.
    x_axis_label : str
        The label of the x-axis. Options are:
            - m/z, Th
            - Inversed IM, V·s·cm\u207B\u00B2
    y_axis_label : str
        The label of the y-axis. Options are:
            - m/z, Th
            - Inversed IM, V·s·cm\u207B\u00B2
    z_axis_label : str
        Should not be set for a 2D scatterplot / heatmap. Default is "Intensity".
    title : str
        The title of the plot. Default is "".
    width : int
        Plot width. Default is 450.
    height : int
        Plot height. Default is 450.
    background_color : str
        The background color of the plot. Default is "black".
    precursor_color : str
        The color of the precursor marker. Default is "green".
    precursor_size : int
        The size of the precursor marker. Default is 15.
    colormap : str
        The colormap name. Default is "fire".
    **kwargs
        Additional keyword arguments to be passed to customization of the hv.Scatter object.

    Returns
    -------
    hv.Points
        A scatter plot projected on the 2 dimensions with markered position of the precursor.

    """
    labels = {
        'm/z, Th': "mz_values",
        'Inversed IM, V·s·cm\u207B\u00B2': "mobility_values",
        'Intensity': "intensity_values",
    }
    x_dimension = labels[x_axis_label]
    y_dimension = labels[y_axis_label]
    z_dimension = labels[z_axis_label]

    opts_ms1 = dict(
        width=width,
        height=height,
        title=title,
        xlabel=x_axis_label,
        ylabel=y_axis_label,
        bgcolor=background_color,
        hooks=[_change_plot],
        **kwargs
    )
    dmap = hv.DynamicMap(
        hv.Points(
            df,
            [x_dimension, y_dimension],
            z_dimension
        )
    )

    agg = rasterize(
        dmap,
        width=width,
        height=height,
        aggregator='sum'
    )
    fig = dynspread(
        shade(
            agg,
            cmap=colormap
        )
    ).opts(plot=opts_ms1)

    if mz and im:
        if x_dimension == 'mz_values' and y_dimension == 'mobility_values':
            precursor = hv.Points((mz, im)).opts(
                marker='x',
                size=precursor_size,
                color=precursor_color,
            )
        else:
            precursor = hv.Points((im, mz)).opts(
                marker='x',
                size=precursor_size,
                color=precursor_color,
            )
        return fig * precursor
    return fig


def _change_plot(plot, element):
    plot.state.toolbar.logo = None


def plot_line(
    timstof_data,  # alphatims.bruker.TimsTOF object
    selected_indices: np.ndarray,
    x_axis_label: str,
    colorscale_qualitative: str,
    title: str = "",
    y_axis_label: str = "intensity",
    remove_zeros: bool = False,
    trim: bool = True,
    height: int = 400
) -> go.Figure:
    """Plot an XIC, mobilogram or spectrum as a lineplot.

    Parameters
    ----------
    timstof_data : alphatims.bruker.TimsTOF object
        An alphatims.bruker.TimsTOF data object.
    selected_indices : np.ndarray
        The raw indices that are selected for this plot. These are typically obtained by slicing the TimsTOF data object with e.g. data[..., "raw"].
    x_axis_label : str
        The label of the x-axis. Options are:
            - mz
            - rt
            - mobility
    y_axis_label : str
        Should not be set for a 1D line plot. Default is "intensity".
    title : str
        The title of the plot. Default is "".
    remove_zeros : bool
        If True, zeros are removed. Note that a line plot connects consecutive points, which can lead to misleading plots if non-zeros are removed. If False, use the full range of the appropriate dimension of the timstof_data. Default is False.
    trim : bool
        If True, zeros on the left and right are trimmed. Default is True.
    height : int
        Plot height. Default is 400.

    Returns
    -------
    plotly.graph_objects.Figure object
        A lne plot showing an XIC, mobilogram or spectrum.

    """
    axis_dict = {
        "mz": "m/z, Th",
        "rt": "RT, min",
        "mobility": "Inversed IM, V·s·cm\u207B\u00B2",
        "intensity": "Intensity",
    }
    x_axis_label = axis_dict[x_axis_label]
    y_axis_label = axis_dict[y_axis_label]
    labels = {
        'm/z, Th': "mz_values",
        'RT, min': "rt_values",
        'Inversed IM, V·s·cm\u207B\u00B2': "mobility_values",
    }
    x_dimension = labels[x_axis_label]
    intensities = timstof_data.bin_intensities(selected_indices, [x_dimension])

    if x_dimension == "mz_values":
        x_ticks = timstof_data.mz_values
        plot_title = "Spectrum"
    elif x_dimension == "mobility_values":
        x_ticks = timstof_data.mobility_values
        plot_title = "Mobilogram"
    elif x_dimension == "rt_values":
        x_ticks = timstof_data.rt_values / 60
        plot_title = "XIC"
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
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x_ticks,
            y=intensities,
            mode='lines',
            text=[f'{x_axis_label}'.format(i + 1) for i in range(len(x_ticks))],
            hovertemplate='<b>%{text}:</b> %{x};<br><b>Intensity:</b> %{y}.',
            name=" ",
            marker=dict(color=getattr(px.colors.qualitative, colorscale_qualitative)[0])
        )
    )

    fig.update_layout(
        title=dict(
            text=plot_title,
            font=dict(
                size=16,
            ),
            x=0.5,
            xanchor='center',
            yanchor='top'
        ),
        xaxis=dict(
            title=x_axis_label,
            titlefont_size=14,
            tickmode='auto',
            tickfont_size=14,
        ),
        yaxis=dict(
            title=y_axis_label,
        ),
        template="plotly_white",
        height=height,
        hovermode="x"
    )

    return fig


def plot_mass_spectra(
    data: pd.DataFrame,
    title: str,
    predicted: tuple = (),
    template: str = "plotly_white",
    spectrum_color: str = 'grey',
    b_ion_color: str = 'red',
    y_ion_color: str = 'blue',
    spectrum_line_width: float = 1.5,
    height: int = 520,
) -> go.Figure:
    """Plot the mass spectrum with a mass error plot for each ion and annotated peptide sequence as subplots.

    Parameters
    ----------
    data : pd.DataFrame
        The dataframe containing spectrum information such as 'mz_values', 'intensity_values', 'ions'.
    title : str
        The title of the plot.
    predicted : tuple
        The tuple containing values of the predicted FragmentMz, RelativeIntensity and ions in the form of:
        (predicted_df.FragmentMz, predicted_df.RelativeIntensity, predicted_df.ions). Default: empty tuple.
    spectrum_color : str
        The color of the mass spectrum. Default is 'grey'.
    b_ion_color : str
        The color of the b-ions. Default is 'red'.
    y_ion_color : str
        The color of the y-ions. Default is 'blue'.
    spectrum_line_width: float
        The width of the spectrum peaks. Default is 1.5.

    Returns
    -------
    plotly.graph_objects.Figure object
        The ms2 spectum plot with the mass error plot for each ion and the annotated peptide sequence as subplots.

    """
    fig = go.Figure()

    if predicted:
        import sklearn.preprocessing
        scaled_int = sklearn.preprocessing.MinMaxScaler((0, 100)).fit_transform(data.intensity_values.values.reshape(-1, 1))
        data['intensity_values'] = scaled_int.reshape(1, -1)[0]

    fig.add_trace(
        go.Scatter(
            x=data[data.ions == '-'].mz_values,
            y=data[data.ions == '-'].intensity_values,
            mode='markers',
            marker=dict(color=spectrum_color, size=1),
            hovertext=data[data.ions == '-'].ions,
            hovertemplate='<b>m/z:</b> %{x};<br><b>Intensity:</b> %{y};<br><b>Ion:</b> %{hovertext}.',
            name='',
            showlegend=False
        )
    )
    # y-ions
    data_y_ions = data[data.ions.str.contains('y')]
    fig.add_trace(
        go.Scatter(
            x=data_y_ions.mz_values,
            y=data_y_ions.intensity_values,
            mode='markers',
            marker=dict(color=y_ion_color, size=1),
            hovertext=data_y_ions.ions,
            hovertemplate='<b>m/z:</b> %{x};<br><b>Intensity:</b> %{y};<br><b>Ion:</b> %{hovertext}.',
            name='',
            showlegend=False
        )
    )
    # b-ions
    data_b_ions = data[data.ions.str.contains('b')]
    fig.add_trace(
        go.Scatter(
            x=data_b_ions.mz_values,
            y=data_b_ions.intensity_values,
            mode='markers',
            marker=dict(color=b_ion_color, size=1),
            hovertext=data_b_ions.ions,
            hovertemplate='<b>m/z:</b> %{x};<br><b>Intensity:</b> %{y};<br><b>Ion:</b> %{hovertext}.',
            name='',
            showlegend=False
        )
    )

    if predicted:
        predicted_b_ions_ind = predicted[2][predicted[2].str.contains('b')].index
        fig.add_trace(
            go.Scatter(
                x=predicted[0][predicted_b_ions_ind],
                y=predicted[1][predicted_b_ions_ind],
                mode='markers',
                opacity=0.7,
                marker=dict(color=b_ion_color, size=1),
                hovertext=predicted[2][predicted_b_ions_ind],
                hovertemplate='<b>m/z:</b> %{x};<br><b>Intensity:</b> %{y};<br><b>Ion:</b> %{hovertext}.',
                name='',
                showlegend=False
            )
        )
        predicted_y_ions_ind = list(set(predicted[2].index).difference(predicted_b_ions_ind))
        fig.add_trace(
            go.Scatter(
                x=predicted[0][predicted_y_ions_ind],
                y=predicted[1][predicted_y_ions_ind],
                mode='markers',
                opacity=0.7,
                marker=dict(color=y_ion_color, size=1),
                hovertext=predicted[2][predicted_y_ions_ind],
                hovertemplate='<b>m/z:</b> %{x};<br><b>Intensity:</b> %{y};<br><b>Ion:</b> %{hovertext}.',
                name='',
                showlegend=False
            )
        )
        fig.update_layout(
            shapes=[
                dict(
                    type='line',
                    xref='x',
                    yref='y',
                    x0=predicted[0][i],
                    y0=0,
                    opacity=0.7,
                    x1=predicted[0][i],
                    y1=predicted[1][i],
                    line=dict(
                        color=b_ion_color if 'b' in predicted[2][i] else y_ion_color,
                        width=spectrum_line_width
                    )
                ) for i in range(len(predicted[0]))
            ],
        )

    fig.update_layout(
        shapes=[
            dict(
                type='line',
                xref='x',
                yref='y',
                x0=data.loc[i, 'mz_values'],
                y0=0,
                x1=data.loc[i, 'mz_values'],
                y1=data.loc[i, 'intensity_values'],
                line=dict(
                    color=b_ion_color if 'b' in data.loc[i, 'ions'] else (y_ion_color if 'y' in data.loc[i, 'ions'] else spectrum_color),
                    width=spectrum_line_width
                )
            ) for i in data.index
        ]
    )

    fig.update_layout(
        template=template,
        xaxis=dict(
            visible=True,
            title='m/z, Th',
        ),
        legend=dict(
            orientation="h",
            x=1,
            xanchor="right",
            yanchor="bottom",
            y=1.01
        ),
        hovermode="x",
        height=height,
        title=dict(
            text=title,
            yanchor='bottom'
        ),
        yaxis=dict(title='Intensity') if not predicted else dict(
            title='Relative intensity, %',
            ticktext=["100", "50", "0", "50", "100"],
            tickvals=[-100, -50, 0, 50, 100],
        ),
    )
    return fig


def plot_complex_ms_plot(
    data: pd.DataFrame,
    title: str,
    sequence: str,
    predicted: tuple = (),
    spectrum_color: str = 'grey',
    template: str = "plotly_white",
    b_ion_color: str = 'red',
    y_ion_color: str = 'blue',
    spectrum_line_width: float = 1.5,
    font_size_seq: int = 14,
    font_size_ion: int = 10,
    height: int = 520
) -> go.Figure:
    """Plot the mass spectrum with a mass error plot for each ion and annotated peptide sequence as subplots.

    Parameters
    ----------
    data : pd.DataFrame
        The dataframe containing spectrum information such as 'mz_values', 'intensity_values', 'ions'.
    title : str
        The title of the plot.
    sequence: str
        The peptide sequence.
    predicted : tuple
        The tuple containing values of the predicted FragmentMz, RelativeIntensity and ions in the form of:
        (predicted_df.FragmentMz, predicted_df.RelativeIntensity, predicted_df.ions). Default: empty tuple.
    spectrum_color : str
        The color of the mass spectrum. Default is 'grey'.
    b_ion_color : str
        The color of the b-ions. Default is 'red'.
    y_ion_color : str
        The color of the y-ions. Default is 'blue'.
    spectrum_line_width: float
        The width of the spectrum peaks. Default is 1.5.
    font_size_seq: int
        The font size of the peptide sequence letters. Default is 14.
    font_size_ion: int
        The font size of the ion letters. Default is 10.
    height: int
        The height of the plot. Default is 520.

    Returns
    -------
    plotly.graph_objects.Figure object
        The ms2 spectum plot with the mass error plot for each ion and the annotated peptide sequence as subplots.

    """
    fig = plot_mass_spectra(data=data, title=title, predicted=predicted, spectrum_color=spectrum_color, template=template, b_ion_color=b_ion_color, y_ion_color=y_ion_color, spectrum_line_width=spectrum_line_width, height=height)

    fig_common = plotly.subplots.make_subplots(
        rows=6, cols=3, shared_xaxes=True,
        figure=fig,
        specs=[
          [{"rowspan": 3, "colspan": 3}, None, None],
          [None, None, None],
          [None, None, None],
          [{"rowspan": 2, "colspan": 3}, None, None],
          [None, None, None],
          [{"colspan": 3}, None, None]
        ],
        vertical_spacing=0.09,

    )

    # add a second plot
    data_b_ions = data[data.ions.str.contains('b')]
    fig_common.add_trace(
        go.Scatter(
            x=data_b_ions.mz_values,
            y=data_b_ions.mass_dev_ppm.round(4),
            hovertext=data_b_ions.ions,
            hovertemplate='<b>m/z:</b> %{x};<br><b>Mass error(ppm):</b> %{y};<br><b>Ion:</b> %{hovertext}.',
            mode='markers',
            marker=dict(
                color=b_ion_color,
                size=6
            ),
            name='b ions'
        ),
        row=4, col=1
    )
    data_y_ions = data[data.ions.str.contains('y')]
    fig_common.add_trace(
        go.Scatter(
            x=data_y_ions.mz_values,
            y=data_y_ions.mass_dev_ppm.round(4),
            hovertext=data_y_ions.ions,
            hovertemplate='<b>m/z:</b> %{x};<br><b>Mass error(ppm):</b> %{y};<br><b>Ion:</b> %{hovertext}.',
            mode='markers',
            marker=dict(
                color=y_ion_color,
                size=6
            ),
            name='y ions'
        ),
        row=4, col=1
    )

    fig_common.update_yaxes(title_text="Error, ppm", row=4, col=1)
    fig_common.update_xaxes(range=[data.mz_values.min()-10, data.mz_values.max()+10], row=4, col=1)
    fig_common.update_xaxes(range=[data.mz_values.min()-10, data.mz_values.max()+10], row=1, col=1)

    bions = alphaviz.preprocessing.get_identified_ions(data.ions, sequence, 'b')
    yions = alphaviz.preprocessing.get_identified_ions(data.ions, sequence, 'y')

    sl = len(sequence)
    distance_from_side = (data.mz_values.max() - data.mz_values.min()) * 2/8
    distance = np.linspace(data.mz_values.min()+distance_from_side, data.mz_values.max()-distance_from_side, sl+1)
    for i, aa in enumerate(sequence):
        fig_common.add_annotation(
            dict(
                text=aa,
                x=distance[i],
                y=0,
                showarrow=False,
                font_size=font_size_seq,
                yshift=1, align='center'
            ),
            row=6,
            col=1
        )
    for i, b in enumerate(bions):
        if b:
            fig_common.add_trace(
                go.Scatter(
                    x=[distance[i], distance[i] + (distance[i+1] - distance[i])/2, distance[i] + (distance[i+1] - distance[i])/2],
                    y=[0.7, 0.7, 0],
                    mode="lines",
                    showlegend=False,
                    marker_color=b_ion_color,
                    line_width=spectrum_line_width,
                    hoverinfo='skip'
                ),
                row=6,
                col=1
            )
            fig_common.add_annotation(
                dict(
                    text="b{}".format(str(i+1)),
                    x=distance[i] + (distance[i+1] - distance[i])/4,
                    y=1.1,
                    showarrow=False,
                    font_size=font_size_ion
                ),
                row=6,
                col=1
            )
    for i, y in enumerate(yions):
        if y:
            fig_common.add_trace(
                go.Scatter(
                    x=[distance[i], distance[i] - (distance[i+1] - distance[i])/2, distance[i] - (distance[i+1] - distance[i])/2],
                    y=[-0.7, -0.7, 0],
                    mode="lines",
                    showlegend=False,
                    marker_color=y_ion_color,
                    line_width=spectrum_line_width,
                    hoverinfo='skip'
                ),
                row=6,
                col=1
            )
            fig_common.add_annotation(
                dict(
                    text="y{}".format(str(sl-i)),
                    x=distance[i] - (distance[i+1] - distance[i])/4,
                    y=-1.1,
                    showarrow=False,
                    font_size=font_size_ion
                ),
                row=6,
                col=1
            )
    fig_common.update_yaxes(
        visible=False,
        range=(-1.1, 1.1),
        row=6,
        col=1
    )
    fig_common.update_xaxes(
        visible=False,
        row=6,
        col=1
    )
    fig_common.update_xaxes(matches='x')
    return fig_common


def plot_mass_error(
    df: pd.DataFrame,
    x_axis_label: str,
    y_axis_label: str,
    plot_title: str,
    mz_tol: float = None
) -> go.Figure:
    """Create a density plot superimposed on the scatter plot together with the 1D distributions of both variables as marginal histograms.

    Parameters
    ----------
    df : pd.DataFrame
        The data frame containing the data.
    x_axis_label : str
        The label of the x-axis.
    y_axis_label : str
        The label of the y-axis.
    plot_title : str
        The title of the plot.

    Returns
    -------
    plotly.graph_objects.Figure object
        Superimposed density and scatter plots with the values distribution of both axes as marginal histograms.

    """
    fig = go.Figure()
    fig.add_trace(
        go.Histogram2dContour(
            x=df[x_axis_label].values,
            y=df[y_axis_label].values,
            colorscale='Blues',
            ncontours=5,
            name=" ",
            contours=dict(
                showlabels=False,
                coloring='fill'
            ),
            hoverinfo='none',
        )
    )
    fig.add_trace(
        go.Scatter(
            x=df[x_axis_label].values,
            y=df[y_axis_label].values,
            mode='markers',
            marker=dict(
                color='rgba(0,0,0,0.3)',
                size=3,
                opacity=0.2
            ),
            name=" ",
            # hoverinfo='none',
        )
    )
    if mz_tol:
        fig.add_trace(
            go.Scatter(
                x=[df[x_axis_label].values.min(), df[x_axis_label].values.max()],
                y=[mz_tol, mz_tol],
                mode='lines',
                line=dict(
                    color='darkred',
                    width=2,
                    dash='dash'
                ),
                showlegend=False,
                hoverinfo='none',
            )
        )
        fig.add_trace(
            go.Scatter(
                x=[df[x_axis_label].values.min(), df[x_axis_label].values.max()],
                y=[-mz_tol, -mz_tol],
                mode='lines',
                line=dict(
                    color='darkred',
                    width=2,
                    dash='dash'
                ),
                showlegend=False,
                hoverinfo='none',
            )
        )

    fig.update_layout(
        autosize=False,
        title=dict(
            text=plot_title,
            font=dict(
                size=16,
            ),
            x=0.5,
            xanchor='center',
            yanchor='top'
        ),
        height=500,
        width=500,
        bargap=0,
        hovermode='closest',
        showlegend=False,
        template='plotly_white'
    )
    return fig


def plot_peptide_distr(
    df: pd.DataFrame,
    x_axis_label: str,
    plot_title: str
) -> go.Figure:
    """Create a distribution plot in conjuction with boxplot.

    Parameters
    ----------
    df : pd.DataFrame
        The data frame containing the data.
    x_axis_label : str
        The label of the x-axis.
    plot_title : str
        The title of the plot.

    Returns
    -------
    plotly.graph_objects.Figure object
        A distribution plot represented as a histogram with a boxplot.

    """
    fig = go.Figure()

    fig.add_trace(
        go.Histogram(
            x=df[x_axis_label],
            xaxis='x',
            nbinsy=50,
            marker=dict(
                color='rgb(198,219,239)'
            ),
            showlegend=False,
        )
    )

    fig.add_trace(
        go.Box(
            x=df[x_axis_label],
            yaxis='y2',
            marker_color='rgb(198,219,239)',
            name='',
            showlegend=False,
        )
    )

    fig.update_layout(
        title=dict(
            text=plot_title,
            font=dict(
                size=16,
            ),
            x=0.5,
            xanchor='center',
            yanchor='top'
        ),
        xaxis=dict(
            zeroline=True,
            domain=[0, 0.85],
            showgrid=True,
            title=x_axis_label,

        ),
        yaxis=dict(
            zeroline=True,
            domain=[0, 0.85],
            showgrid=True,
            title='Count',
        ),
        yaxis2=dict(
            zeroline=True,
            domain=[0.85, 1],
            showgrid=True,
            title=''
        ),
        bargroupgap=0.1,
        template='plotly_white',
        showlegend=False,
        height=400,
        width=600,
    )
    return fig


def plot_elution_heatmap(
    df: pd.DataFrame,
    title: str = "",
    width: int = 250,
    height: int = 250,
    background_color: str = "black",
    colormap: str = "fire",
    **kwargs
):
    """Create a heatmap showing a correlation of retention time  and ion mobility with color coding for signal intensity.

    Parameters
    ----------
    df : pandas Dataframe
        A dataframe obtained by slicing an alphatims.bruker.TimsTOF object.
    title: str
        The title of the plot. Default: "".
    width : int
        The width of the plot. Default: 250.
    height : int
        The height of the plot. Default: 250.
    background_color : str
        The background color of the plot. Default: "black".
    colormap : str
        The name of the colormap in Plotly. Default: "fire".

    Returns
    -------
    a Plotly scatter plot
        The scatter plot showing the correlation of retention time  and ion mobility with color coding for signal intensity.
    """
    labels = {
        'RT, min': "rt_values",
        'Inversed IM, V·s·cm\u207B\u00B2': "mobility_values",
        'Intensity': "intensity_values",
    }
    x_axis_label = "RT, min"
    y_axis_label = "Inversed IM, V·s·cm\u207B\u00B2"
    z_axis_label = "Intensity"

    x_dimension = labels[x_axis_label]
    y_dimension = labels[y_axis_label]
    z_dimension = labels[z_axis_label]

    df["rt_values"] /= 60

    opts_ms1 = dict(
        width=width,
        height=height,
        title=title,
        xlabel=x_axis_label,
        ylabel=y_axis_label,
        bgcolor=background_color,
        framewise=True,
        axiswise=True,
        fontsize={'title': 12},
        **kwargs
    )
    dmap = hv.DynamicMap(
        hv.Points(
            df,
            [x_dimension, y_dimension],
            z_dimension
        )
    )
    agg = rasterize(
        dmap,
        width=width,
        height=height,
        aggregator='sum'
    )
    fig = dynspread(
        shade(
            agg,
            cmap=colormap
        )
    ).opts(plot=opts_ms1)

    return fig


def plot_elution_profile_heatmap(
    timstof_data,
    peptide_info: dict,
    mass_dict: dict,
    mz_tol: int = 50,
    rt_tol: int = 30,
    im_tol: int = 0.05,
    title: str = "",
    n_cols: int = 5,
    # width: int = 180,
    height: int = 400,
    background_color: str = "black",
    colormap: str = 'fire',
    **kwargs
):
    """Plot an elution profile for the specified precursor and all his identified fragments as heatmaps in the
    retention time/ion mobility dimensions.

    Parameters
    ----------
    timstof_data : alphatims.bruker.TimsTOF
        An alphatims.bruker.TimsTOF data object.
    peptide_info : dict
        Peptide information including sequence, fragments' patterns, rt, mz and im values.
    mass_dict : dict
        The basic mass dictionaty with the masses of all amino acids and modifications.
    mz_tol: float
        The mz tolerance value. Default: 50 ppm.
    rt_tol: float
        The rt tolerance value. Default: 30 ppm.
    im_tol: float
        The im tolerance value. Default: 0.05 ppm.
    title : str
        The title of the plot. Default: "".
    n_cols: int
        The number of the heatmaps plotted per row. Default: 5.
    # width : int
    #     The width of the plot. Default: 180.
    height : int
        The height of the plot. Default: 180.

    Returns
    -------
    a Bokeh heatmap plots
        The elution profile heatmap plots in retention time and ion mobility dimensions
        for the specified peptide and all his fragments.
    """
    # predict the theoretical fragments using the Alphapept get_fragmass() function.
    frag_masses, frag_type = alphaviz.utils.get_fragmass(
        parsed_pep=alphaviz.utils.parse(peptide_info['sequence']),
        mass_dict=mass_dict
    )
    peptide_info['fragments'] = {
        (f"b{key}" if key > 0 else f"y{-key}"): value for key, value in zip(frag_type, frag_masses)
    }

    # slice the data using the rt_tol, im_tol and mz_tol values
    rt_slice = slice(peptide_info['rt'] - rt_tol, peptide_info['rt'] + rt_tol)
    im_slice = slice(peptide_info['im'] - im_tol, peptide_info['im'] + im_tol)
    prec_mz_slice = slice(peptide_info['mz'] / (1 + mz_tol / 10**6), peptide_info['mz'] * (1 + mz_tol / 10**6))

    # create an elution profile for the precursor
    precursor_indices = timstof_data[
        rt_slice,
        im_slice,
        0,
        prec_mz_slice,
        'raw'
    ]

    common_plot = plot_elution_heatmap(
        timstof_data.as_dataframe(precursor_indices),
        title="precursor",
        # width=width,
        height=height,
        background_color=background_color,
        colormap=colormap,
        # ylim=(im_slice.start, im_slice.stop),
        # y_axis_label="RT, min",
        # x_axis_label="Inversed IM, V·s·cm\u207B\u00B2",
        **kwargs
    )

    # create elution profiles for all fragments
    for frag, frag_mz in peptide_info['fragments'].items():
        fragment_data_indices = timstof_data[
            rt_slice,
            im_slice,
            prec_mz_slice,
            slice(frag_mz / (1 + mz_tol / 10**6), frag_mz * (1 + mz_tol / 10**6)),
            'raw'
        ]
        if len(fragment_data_indices) > 0:
            common_plot += plot_elution_heatmap(
                timstof_data.as_dataframe(fragment_data_indices),
                title=f"{frag} ({round(frag_mz, 3)})",
                # width=width,
                height=height,
                background_color=background_color,
                colormap=colormap,
                # y_axis_label="RT, min",
                # x_axis_label="Inversed IM, V·s·cm\u207B\u00B2",
                # ylim=(im_slice.start, im_slice.stop),
                **kwargs
            )
    try:
        return common_plot.cols(n_cols)
    except AttibuteError:
        return common_plot


def plot_elution_line(
    timstof_data,
    selected_indices: np.ndarray,
    label: str,
    marker_color: str,
    remove_zeros: bool = False,
    trim: bool = True,
):
    """Plot an XIC as a lineplot.

    Parameters
    ----------
    timstof_data : alphatims.bruker.TimsTOF
        An alphatims.bruker.TimsTOF data object.
    selected_indices : np.ndarray
        The raw indices that are selected for this plot.
    label : str
        The label for the line plot.
    remove_zeros : bool
        If True, zeros are removed. Default: False.
    trim : bool
        If True, zeros on the left and right are trimmed. Default: True.

    Returns
    -------
    a Plotly line plot
        The XIC line plot.
    """
    axis_dict = {
        "rt": "RT, min",
        "intensity": "Intensity",
    }
    x_axis_label = axis_dict["rt"]
    # y_axis_label = axis_dict["intensity"]
    labels = {
        'RT, min': "rt_values",
    }
    x_dimension = labels[x_axis_label]
    intensities = timstof_data.bin_intensities(selected_indices, [x_dimension])
    x_ticks = timstof_data.rt_values / 60

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
        hovertemplate='<b>%{text}:</b> %{x};<br><b>Intensity:</b> %{y}.',
        name=label,
        marker=marker_color,
    )
    return trace


def plot_elution_profile(
    raw_data,
    peptide_info: dict,
    mass_dict: dict,
    colorscale_qualitative: str,
    colorscale_sequential: str,
    calculate_fragment_masses: bool = True,
    mz_tol: float = 50,
    rt_tol: float = 30,
    im_tol: float = 0.05,
    title: str = "",
    # width: int = 900,
    height: int = 400,
):
    """Plot an elution profile plot for the specified precursor and all
    his identified fragments.

    Parameters
    ----------
    raw_data : alphatims.bruker.TimsTOF
        An alphatims.bruker.TimsTOF data object.
    peptide_info : dict
        Peptide information including sequence, fragments' patterns, rt, mz
        and im values.
    mass_dict : dict
        The basic mass dictionaty with the masses of all amino acids
        and modifications.
    calculate_fragment_masses : bool
        Whether the fragment masses should be calculated inside the function or
        they are already provided in the peptide_info dict. Default: True.
    mz_tol: float
        The mz tolerance value. Default: 50 ppm.
    rt_tol: float
        The rt tolerance value. Default: 30 ppm.
    im_tol: float
        The im tolerance value. Default: 0.05 ppm.
    title : str
        The title of the plot.
    # width : int
    #     The width of the plot. Default: 900.
    height : int
        The height of the plot. Default: 400.

    Returns
    -------
    a Plotly line plot
        The elution profile plot in retention time dimension for the specified peptide and all his fragments.
    """
    import alphaviz.utils

    x_axis_label = "rt"
    y_axis_label = "intensity"

    if calculate_fragment_masses:
        # predict the theoretical fragments using the Alphapept get_fragmass() function.
        frag_masses, frag_type = alphaviz.utils.get_fragmass(
            parsed_pep=alphaviz.utils.parse(peptide_info['sequence']),
            mass_dict=mass_dict
        )
        peptide_info['fragments'] = {
            (f"b{key}" if key > 0 else f"y{-key}"): value for key, value in zip(frag_type, frag_masses)
        }

    # slice the data using the rt_tol, im_tol and mz_tol values
    rt_slice = slice(peptide_info['rt'] - rt_tol, peptide_info['rt'] + rt_tol)
    im_slice = slice(peptide_info['im'] - im_tol, peptide_info['im'] + im_tol)
    prec_mz_slice = slice(peptide_info['mz'] / (1 + mz_tol / 10**6), peptide_info['mz'] * (1 + mz_tol / 10**6))

    if len(peptide_info['fragments'].values()) + 1 <= len(getattr(px.colors.qualitative, colorscale_qualitative)):
        colors_set = getattr(px.colors.qualitative, colorscale_qualitative)
    else:
        colors_set = px.colors.sample_colorscale(colorscale_sequential, samplepoints=len(peptide_info['fragments'].values()) + 1)

    # create an elution profile for the precursor
    precursor_indices = raw_data[
        rt_slice,
        im_slice,
        0,
        prec_mz_slice,
        'raw'
    ]
    fig = go.Figure()
    fig.add_trace(
        plot_elution_line(
            raw_data,
            precursor_indices,
            remove_zeros=True,
            # label=f"precursor ({round(peptide_info['mz'], 3)})",
            label='precursor',
            marker_color=dict(color=colors_set[0])
        )
    )
    # create elution profiles for all fragments
    for ind, (frag, frag_mz) in enumerate(peptide_info['fragments'].items()):
        fragment_data_indices = raw_data[
            rt_slice,
            im_slice,
            prec_mz_slice,
            slice(frag_mz / (1 + mz_tol / 10**6), frag_mz * (1 + mz_tol / 10**6)),
            'raw'
        ]
        if len(fragment_data_indices) > 0:
            fig.add_trace(
                plot_elution_line(
                    raw_data,
                    fragment_data_indices,
                    remove_zeros=True,
                    label=f"{frag} ({round(frag_mz, 3)})",
                    marker_color=dict(color=colors_set[ind+1])
                )
            )

    fig.update_layout(
        title=dict(
            text=title,
            font=dict(
                size=16,
            ),
            x=0.5,
            xanchor='center',
            yanchor='top'
        ),
        xaxis=dict(
            title=x_axis_label,
            titlefont_size=14,
            tickmode='auto',
            tickfont_size=14,
        ),
        yaxis=dict(
            title=y_axis_label
        ),
        legend=dict(
            orientation='h',
            yanchor='bottom',
            y=-1,
            xanchor='right',
            x=0.95,
            font=dict(
                # family="Courier",
                size=11,
                color='black'
            ),
        ),
        template='plotly_white',
        # width=width,
        height=height,
        hovermode='x unified',
        showlegend=True
    )
    return fig


def plot_pept_per_protein_barplot(
    df: pd.DataFrame,
    x_axis_label: str,
    plot_title: str
) -> go.Figure:
    """Create a barplot for the number of peptides identified per protein.

    Parameters
    ----------
    df : pd.DataFrame
        The data frame containing the data.
    x_axis_label : str
        The label of the x-axis containing information about the number of identified peptides per protein.
    plot_title : str
        The title of the plot.

    Returns
    -------
    plotly.graph_objects.Figure object
        A distribution barplot.

    """
    df = df.copy()
    df['pept_per_prot'] = df[x_axis_label].apply(lambda x: str(x) if x < 5 else '>5')

    fig = go.Figure()

    fig.add_trace(
        go.Bar(
            x=df.pept_per_prot.value_counts().sort_index().index,
            y=df.pept_per_prot.value_counts().sort_index().values,
            marker=dict(
                color='rgb(198,219,239)'
            ),
            text=[f'{each:.2f}' for each in df.pept_per_prot.value_counts(normalize=True).sort_index().values],
            textfont={
                'size': 8,
                'color': 'green'
            },
        )
    )

    fig.update_layout(
        title=dict(
            text=plot_title,
            font=dict(
                size=16,
            ),
            x=0.5,
            xanchor='center',
            yanchor='top'
        ),
        xaxis=dict(
            zeroline=True,
            showgrid=True,
            title='Number of peptide',

        ),
        yaxis=dict(
            zeroline=True,
            showgrid=True,
            title='Count',
        ),
        template='plotly_white',
        height=400,
        width=400,
    )
    return fig


def export_svg(obj, filename='test', width=500, height=500):
    plot_state = hv.renderer('bokeh').get_plot(obj).state
    plot_state.output_backend = 'svg'
    export_svgs(plot_state, filename=filename, width=width, height=height)
