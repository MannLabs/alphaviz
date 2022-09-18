import os
import torch
import plotly

import pandas as pd
import numpy as np

import plotly.graph_objects as go

from .msplot_utils import _plot_scatter

color_map:dict = {
    '-': 'lightgrey', # '-' means umnatched
    'b': 'blue', 'y': 'red',
}

class MS2_Plot:
    vertical_spacing = 0.05
    template = 'plotly_white'
    plot_height = 600
    def __init__(self, 
        peak_plot_rows = 4,
        mass_err_plot_rows = 1,
        frag_coverage_plot_rows = 1,
        frag_coverage_at_top = True
    ):
        specs = []
        if frag_coverage_at_top:
            specs.append(
                [{
                    "rowspan": frag_coverage_plot_rows, 
                    "colspan": 3
                }, None, None]
            )
            specs.extend(
                [[None, None, None]]*
                (frag_coverage_plot_rows-1)
            )
        specs.append(
            [{
                "rowspan": peak_plot_rows, 
                "colspan": 3
            }, None, None]
        )
        specs.extend(
            [[None, None, None]]*
            (peak_plot_rows-1)
        )
        specs.append(
            [{
                "rowspan": mass_err_plot_rows, 
                "colspan": 3
            }, None, None]
        )
        specs.extend(
            [[None, None, None]]*
            (mass_err_plot_rows-1)
        )
        if not frag_coverage_at_top:
            specs.append(
                [{
                    "rowspan": frag_coverage_plot_rows, 
                    "colspan": 3
                }, None, None]
            )
            specs.extend(
                [[None, None, None]]*
                (frag_coverage_plot_rows-1)
            )

        if frag_coverage_at_top:
            (
                frag_cov_row,
                peak_row, 
                mass_err_row,
            ) = np.cumsum([
                1, frag_coverage_plot_rows,
                peak_plot_rows,
            ])
        else:
            (
                peak_row, 
                mass_err_row,
                frag_cov_row, 
            ) = np.cumsum([
                1, peak_plot_rows,
                mass_err_plot_rows
            ])

        self.specs = specs
        self.peak_row = peak_row
        self.mass_err_row = mass_err_row
        self.frag_cov_row = frag_cov_row
        self.rows = (
            peak_plot_rows+mass_err_plot_rows+
            frag_coverage_plot_rows
        )

    def plot(self, plot_df, sequence, title, 
        plot_unmatched_peaks=False
    ):
        self._init_plot(title)

        self.peak_plot.plot(
            plot_df, 
            plot_unmatched_peaks=plot_unmatched_peaks
        )
        self.mass_err_plot.plot(
            plot_df,
        )
        self.frag_cov_plot.plot(
            plot_df, sequence
        )
        return self.fig

    def _init_plot(self, title):

        self.fig = plotly.subplots.make_subplots(
            rows=(
                self.rows
            ), cols=3, 
            shared_xaxes=True,
            specs=self.specs,
            vertical_spacing=self.vertical_spacing,
            column_widths=[0.25, 0.5, 0.25]
        )

        self.peak_plot = PeakPlot(
            self.fig, self.peak_row, 
        )
        self.mass_err_plot = MassErrPlot(
            self.fig, self.mass_err_row,
        )
        self.frag_cov_plot = FragCoveragePlot(
            self.fig, self.frag_cov_row
        )

        self.fig.update_layout(
            template=self.template,
            title=dict(
                text=title,
                yanchor='bottom'
            ),
            hovermode='x',
            height=self.plot_height,
        )
        self.fig.update_xaxes(matches='x')
        self.fig.update_yaxes(
            title = 'intensity',
        )

class MassErrPlot:
    def __init__(self, fig_subplots, row):
        """
        See `MS_Plot_Base` Parameters
        """
        self.fig = fig_subplots
        self.row = row
        self.col = 1
        self.hovertemplate = (
            '%{hovertext}<br>'
            '<b>m/z:</b> %{x}<br>'
            '<b>Mass err (ppm):</b> %{y}'
        )

    def plot(self, plot_df):
        self._plot_one_ion(
            plot_df, 'b', color_map['b']
        )
        self._plot_one_ion(
            plot_df, 'y', color_map['y']
        )
        self.fig.update_yaxes(
            title_text='ppm',
            row=self.row,col=self.col
        )
        return self.fig

    def _plot_one_ion(self,
        plot_df, ion_type, color
    ):
        df = plot_df[
            plot_df.ions.str.startswith(ion_type)
        ].query("intensity_values > 0")
        _plot_scatter(
            self.fig,
            df.mz_values, 
            df.mass_dev_ppm.round(4),
            color=color,
            marker_size=5,
            hovertext=df.ions,
            hovertemplate=self.hovertemplate,
            name=f"{ion_type} ions",
            row=self.row,
            col=self.col
        )

class PeakPlot:
    def __init__(self, fig_subplots, row):
        """
        See `MS_Plot_Base` Parameters
        """
        self.fig = fig_subplots
        self.row = row
        self.col = 1
        self.hovertemplate = (
            '<b>%{hovertext}</b><br>'
            '<b>m/z:</b> %{x}<br>'
            '<b>intensity:</b> %{y}'
        )
        self.peak_line_width = 1.5
        
    def plot(self,
        plot_df,
        plot_unmatched_peaks:bool=True,
    )->go.Figure:
        if plot_unmatched_peaks:
            _df = plot_df[plot_df.ions=="-"]
            _plot_scatter(
                self.fig,
                _df.mz_values,
                _df.intensity_values,
                color=color_map['-'],
                marker_size=1,
                hovertext=None,
                hovertemplate=None,
                name='',
                row=self.row,
                col=self.col,
            )
        
        self._plot_one_ion_type_scatter(
            plot_df, 'b'
        )
        self._plot_one_ion_type_scatter(
            plot_df, 'y'
        )

        self._plot_peak_vlines(
            plot_df,
        )

        self._plot_frag_annotations(
            plot_df.query('ions != "-"')
        )
        return self.fig

    def _plot_one_ion_type_scatter(self, 
        plot_df, ion_type
    ):
        df_all = plot_df[
            plot_df.ions.str.startswith(ion_type)
        ]
        _df = df_all.query('intensity_values>0')
        _plot_scatter(
            self.fig,
            _df.mz_values,
            _df.intensity_values,
            color=color_map[ion_type], 
            marker_size=1,
            hovertext=_df.ions,
            hovertemplate=self.hovertemplate,
            name='',
            row=self.row,
            col=self.col,
        )
        _df = df_all.query('intensity_values<0')
        _plot_scatter(
            self.fig,
            _df.mz_values,
            _df.intensity_values,
            color=color_map[ion_type], 
            marker_size=1,
            hovertext=_df.ions,
            hovertemplate=self.hovertemplate,
            name='',
            row=self.row,
            col=self.col,
        )

    def _plot_frag_annotations(self, plot_df):
        df = plot_df.query('intensity_values>0')
        max_inten = df.intensity_values.max()
        yshift=max_inten*0.02
        for mz, inten, ion in df[
            ['mz_values','intensity_values','ions']
        ].values:
            self.fig.add_annotation(
                x=mz, y=inten+yshift,
                text=ion,
                textangle=-90,
                font_size=10,
                row=self.row,
                col=self.col,
            )
        
        neg_ay = max_inten*0.3
        pred_df = plot_df.query('intensity_values<0')
        pred_df = pred_df[~pred_df.ions.isin(set(df.ions))]
        for mz, inten, ion in pred_df[
            ['mz_values','intensity_values','ions']
        ].values:
            self.fig.add_annotation(
                x=mz, y=inten-yshift,
                text=ion,
                textangle=-90,
                font_size=10,
                ay=inten-yshift-neg_ay,
                ayref=f'y{self.row}',
                yref=f'y{self.row}',
                row=self.row,
                col=self.col,
            )

    def _plot_peak_vlines(self,
        plot_df,
    ):
        self.fig.update_layout(
            shapes=[
                dict(
                    type='line',
                    xref=f'x{self.row}',
                    yref=f'y{self.row}',
                    x0=plot_df.loc[i, 'mz_values'],
                    y0=0,
                    x1=plot_df.loc[i, 'mz_values'],
                    y1=plot_df.loc[i, 'intensity_values'],
                    line=dict(
                        color=color_map[plot_df.loc[i,'ions'][0]],
                        width=self.peak_line_width,
                    )
                ) for i in plot_df.index
            ]
        )
        # for i in plot_df.index:
        #     self.fig.add_shape(type='line',
        #         x0=plot_df.loc[i, 'mz_values'],
        #         x1=plot_df.loc[i, 'mz_values'],
        #         y0=0,
        #         y1=plot_df.loc[i, 'intensity_values'],
        #         line_color = color_map[plot_df.loc[i, 'ions'][0]],
        #         line_width = self.peak_line_width,
        #         row=self.row,
        #         col=self.col,
        #     )

class FragCoveragePlot:
    def __init__(self, fig_subplots, row):
        """
        See `MS_Plot_Base` Parameters
        """
        self.fig = fig_subplots
        self.row = row
        self.col = 1
        self.font_size_sequence = 14
        self.font_size_coverage = 8

    def plot(self,
        plot_df,
        sequence,
    ):
        d = (
            plot_df.mz_values.max() - 
            plot_df.mz_values.min()
        ) * 2/8
        aa_x_positions = np.linspace(
            plot_df.mz_values.min()+d, 
            plot_df.mz_values.max()-d, 
            len(sequence)+1
        )

        self._plot_sequence(sequence, aa_x_positions)
        self._plot_coverage_one_frag_type(
            plot_df, sequence, aa_x_positions, 
            'b', color_map['b'],
        )
        self._plot_coverage_one_frag_type(
            plot_df, sequence, aa_x_positions, 
            'y', color_map['y'],
        )
        self.fig.update_yaxes(
            visible=False,
            range=(-1.1, 1.1),
            row=self.row,
            col=self.col,
        )
        self.fig.update_xaxes(
            visible=False,
            row=self.row,
            col=self.col,
        )
        return self.fig

    def _plot_sequence(self,
        sequence, aa_x_positions,
    ):
        for i, aa in enumerate(sequence):
            self.fig.add_annotation(
                dict(
                    text=aa,
                    x=aa_x_positions[i],
                    y=0,
                    showarrow=False,
                    font_size=self.font_size_sequence,
                    yshift=1, align='center'
                ),
                row=self.row,
                col=self.col
            )

    def _plot_coverage_one_frag_type(self,
        plot_df, sequence, 
        aa_x_positions, 
        ion_type,
        color,
    ):
        nAA = len(sequence)
        plot_df = plot_df[
            plot_df.ions.str.startswith(ion_type)
        ].query("intensity_values>0")

        covs = np.zeros(nAA, dtype=np.int64)
        if ion_type in 'abc':
            covs[plot_df.fragment_indices] = 1
        else:
            covs[plot_df.fragment_indices+1] = 1

        def get_positions(ion_type, i):
            if ion_type in 'abc':
                return dict(
                    x=[
                        aa_x_positions[i], 
                        aa_x_positions[i] + (
                            aa_x_positions[i+1] - aa_x_positions[i]
                        )/2, 
                        aa_x_positions[i] + (
                            aa_x_positions[i+1] - aa_x_positions[i]
                        )/2
                    ], 
                    y=[0.6,0.6,0.4],
                    tx=aa_x_positions[i] + (
                        aa_x_positions[i+1] - aa_x_positions[i]
                    )/4,
                    ty=1.0,
                    s=i+1
                )
            else:
                return dict(
                    x=[
                        aa_x_positions[i], 
                        aa_x_positions[i] - (
                            aa_x_positions[i+1] - aa_x_positions[i]
                        )/2, 
                        aa_x_positions[i] - (
                            aa_x_positions[i+1] - aa_x_positions[i]
                        )/2
                    ],
                    y=[-0.6, -0.6, -0.4],
                    tx=aa_x_positions[i] - (
                        aa_x_positions[i+1] - aa_x_positions[i]
                    )/4,
                    ty=-1.0,
                    s=nAA-i
                )
        for i, cov in enumerate(covs):
            if cov:
                pos = get_positions(ion_type, i)
                self.fig.add_trace(
                    go.Scatter(
                        x=pos['x'],
                        y=pos['y'],
                        mode="lines",
                        showlegend=False,
                        marker_color=color,
                        line_width=1,
                        hoverinfo='skip'
                    ),
                    row=self.row,
                    col=self.col,
                )
                self.fig.add_annotation(
                    dict(
                        text=f"{ion_type}{pos['s']}",
                        x=pos['tx'],
                        y=pos['ty'],
                        showarrow=False,
                        font_size=self.font_size_coverage
                    ),
                    row=self.row,
                    col=self.col,
                )
