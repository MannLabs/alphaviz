import os
import torch
import plotly

import pandas as pd
import numpy as np

import plotly.graph_objects as go
import plotly.express as px

from alphatims.bruker import (
    TimsTOF, 
)

from alphaviz.contrib.msplot_utils import (
    _plot_line, _plot_line_fast
)

class XIC_1D_Plot():
    # hovermode = "x" | "y" | "closest" | False | "x unified" | "y unified"
    hovermode = 'closest'
    plot_height = 550
    colorscale_qualitative="Alphabet"
    colorscale_sequential="Viridis"
    theme_template='plotly_white'
    min_frag_mz = 200
    view_dim = 'rt' # or 'im'

    def plot(self,
        tims_data:TimsTOF,
        peptide_info: pd.DataFrame,
        mz_tol: float = 50,
        rt_tol: float = 30,
        im_tol: float = 0.05,
        include_precursor:bool=True,
        include_ms1:bool=True,
    )->go.Figure:
        """Based on `alphaviz.plotting.plot_elution_profile`

        Parameters
        ----------
        peptide_info : pd.DataFrame
            alphaviz peptide_info df

        mz_tol : float, optional
            in ppm, by default 50

        rt_tol : float, optional
            RT tol in seconds, by default 30

        im_tol : float, optional
            mobility tol, by default 0.05

        height : int, optional
            fig height, by default 400

        Returns
        -------
        go.Figure
            plotly Figure object return by 
            `alphaviz.plotting.plot_elution_profile`
        """
        return self.plot_elution_profiles(
            tims_data, peptide_info_list=peptide_info,
            mz_tol=mz_tol, rt_tol=rt_tol, im_tol=im_tol,
            include_precursor=include_precursor,
            include_ms1=include_ms1,
        )

    def _add_trace_df(self,
        fig, row, col, df, view_indices_df, label, marker_color
    ):
        if len(df) == 0: return
        fig.add_trace(
            _plot_line(
                df,
                view_indices_df,
                label=label,
                marker_color=marker_color, 
                view_dim=self.view_dim,
            ),
            row=row, col=col,
        )

    def _add_trace_fast(self,
        fig, row, col, tims_data, selected_indices, view_indices,
        label, marker_color
    ):
        if len(selected_indices) == 0: return
        fig.add_trace(
            _plot_line_fast(
                tims_data, selected_indices,
                view_indices,
                label=label,
                marker_color=marker_color, 
                view_dim=self.view_dim,
            ),
            row=row, col=col,
        )

    def _add_trace(self, 
        fig, row, col, 
        tims_data, selected_indices, view_indices,
        df, view_indices_df, 
        label, marker_color,
    ):
        return self._add_trace_fast(
            fig, row, col, 
            tims_data, selected_indices, view_indices,
            label, marker_color,
        )
        return self._add_trace_df(
            fig, row, col, 
            df, view_indices_df, 
            label, marker_color,
        )

    def _get_view_df_from_indices(self, tims_data, view_indices):
        if self.view_dim == 'rt':
            sliced_view_values = tims_data.rt_values[view_indices]
            return pd.DataFrame({
                'frame_indices': view_indices,
                'rt_values': sliced_view_values/60,
            })
        else:
            sliced_view_values = tims_data.mobility_values[view_indices]
            return pd.DataFrame({
                'scan_indices': view_indices,
                'mobility_values': sliced_view_values/60,
            })

    def _add_elution_profile(self,
        fig:go.Figure,
        row:int, col:int,
        tims_data: TimsTOF,
        peptide_info: pd.DataFrame,
        mz_tol: float = 50,
        rt_tol: float = 30,
        im_tol: float = 0.05,
        include_precursor:bool = True,
        include_ms1:bool = True,
    ):
        """Plot an elution profile plot for the specified precursor and all
        his identified fragments.

        Parameters
        ----------
        raw_data : alphatims.bruker.TimsTOF
            An alphatims.bruker.TimsTOF data object.

        peptide_info : pd.DataFrame
            Peptide information including sequence, fragment patterns, rt,
            and im values.

        mz_tol: float
            The mz tolerance value. Default: 50 ppm.

        rt_tol: float
            The rt tolerance value. Default: 30 ppm.

        im_tol: float
            The im tolerance value. Default: 0.05 ppm.

        Returns
        -------
        a Plotly line plot
            The elution profile plot in retention time dimension for the specified peptide and all his fragments.
        """

        # slice the data using the rt_tol, im_tol and mz_tol values
        rt_slice = slice(
            peptide_info['rt_sec'].values[0] - rt_tol, 
            peptide_info['rt_sec'].values[0] + rt_tol
        )
        im_slice = slice(
            peptide_info['im'].values[0] - im_tol, 
            peptide_info['im'].values[0] + im_tol
        )
        prec_mz_slice = slice(
            peptide_info['precursor_mz'].values[0]
             * (1 - mz_tol / 10**6), 
            peptide_info['precursor_mz'].values[0]
             * (1 + mz_tol / 10**6)
        )

        if len(peptide_info) + 2 <= len(
            getattr(px.colors.qualitative, self.colorscale_qualitative)
        ):
            colors_set = getattr(
                px.colors.qualitative, self.colorscale_qualitative
            )
        else:
            colors_set = px.colors.sample_colorscale(
                self.colorscale_sequential, 
                samplepoints=len(peptide_info) + 2
            )

        ms1_raw_indices = tims_data[
            rt_slice,
            im_slice,
            0,
            :,
            'raw'
        ]

        ms2_raw_indices = tims_data[
            rt_slice,
            im_slice,
            prec_mz_slice,
            :,
            'raw'
        ]
        
        if self.view_dim == 'rt':
            ms1_view_indices = np.sort(np.array(list(set(
                tims_data.convert_from_indices(
                ms1_raw_indices,
                return_frame_indices=True
            )['frame_indices'])), dtype=np.int64))

            ms2_view_indices = np.sort(np.array(list(set(
                tims_data.convert_from_indices(
                ms2_raw_indices,
                return_frame_indices=True
            )['frame_indices'])), dtype=np.int64))
        else:
            ms1_view_indices = np.sort(np.array(list(set(
                tims_data.convert_from_indices(
                ms1_raw_indices,
                return_scan_indices=True
            )['scan_indices'])), dtype=np.int64))

            ms2_view_indices = np.sort(np.array(list(set(
                tims_data.convert_from_indices(
                ms2_raw_indices,
                return_scan_indices=True
            )['scan_indices'])), dtype=np.int64))

        ms1_view_df = self._get_view_df_from_indices(
            tims_data, ms1_view_indices
        )
        ms2_view_df = self._get_view_df_from_indices(
            tims_data, ms2_view_indices
        )

        # create an elution profile for the precursor in ms1
        if include_ms1 and include_precursor:
            ms1_m0_indices = tims_data[
                rt_slice,
                im_slice,
                0,
                prec_mz_slice,
                'raw'
            ]
            self._add_trace_fast(
                fig, row, col,
                tims_data, ms1_m0_indices, ms1_view_indices,
                label=f'MS1 M0 ({peptide_info["precursor_mz"].values[0]:.3f})',
                marker_color=dict(color=colors_set[0])
            )
            # self._add_trace_df(
            #     fig, row, col,
            #     tims_data.as_dataframe(ms1_m0_indices), ms1_view_df,
            #     label=f'MS1 M0 ({peptide_info["precursor_mz"].values[0]:.3f})',
            #     marker_color=dict(color=colors_set[0])
            # )
        # create elution profiles for all fragments
        if include_precursor:
            ms2_m0_indices = tims_data[
                rt_slice,
                im_slice,
                prec_mz_slice,
                prec_mz_slice,
                'raw'
            ]
            self._add_trace_fast(
                fig, row, col,
                tims_data, ms2_m0_indices, ms2_view_indices,
                label=f'MS2 M0 ({peptide_info["precursor_mz"].values[0]:.3f})',
                marker_color=dict(color=colors_set[1])
            )
            # self._add_trace_df(
            #     fig, row, col,
            #     tims_data.as_dataframe(ms2_m0_indices), ms2_view_df,
            #     label=f'MS2 M0 ({peptide_info["precursor_mz"].values[0]:.3f})',
            #     marker_color=dict(color=colors_set[1])
            # )
        for ind, (frag, frag_mz) in enumerate(
            peptide_info[['ion','frag_mz']].values
        ):
            frag_mz = float(frag_mz)
            if frag_mz < self.min_frag_mz: continue
            frag_indices = tims_data[
                rt_slice,
                im_slice,
                prec_mz_slice,
                slice(
                    frag_mz * (1 - mz_tol / 10**6), 
                    frag_mz * (1 + mz_tol / 10**6)
                ),
                'raw'
            ]
            self._add_trace_fast(
                fig, row, col,
                tims_data, frag_indices, ms2_view_indices,
                label=f"{frag} ({frag_mz:.3f})",
                marker_color=dict(color=colors_set[ind+2])
            )
            # self._add_trace_df(
            #     fig, row, col,
            #     tims_data.as_dataframe(frag_indices), ms2_view_df,
            #     label=f"{frag} ({frag_mz:.3f})",
            #     marker_color=dict(color=colors_set[ind+2])
            # )
        
        fig.add_vline(
            peptide_info['rt_sec'].values[0]/60 if 
            self.view_dim == 'rt' else 
            peptide_info['im'].values[0], 
            line_dash="dash", 
            line_color="grey", row=row, col=col,
        )
        return fig

    def plot_elution_profiles(self,
        tims_data:TimsTOF,
        peptide_info_list: pd.DataFrame,
        mz_tol: float = 50,
        rt_tol: float = 30,
        im_tol: float = 0.05,
        include_precursor:bool = True,
        include_ms1:bool = True,
    ):
        if isinstance(peptide_info_list, pd.DataFrame):
            peptide_info_list = [peptide_info_list]

        fig = plotly.subplots.make_subplots(
            rows=len(peptide_info_list), 
            cols=1, 
            shared_xaxes=True,
            x_title='RT (min)' if self.view_dim == 'rt' else 'Mobility',
            y_title='Intensity',
            vertical_spacing=0.2/len(peptide_info_list),
            subplot_titles=[
                _['mod_seq_charge'].values[0] for _ in peptide_info_list
            ]
        )

        for i,peptide_info in enumerate(peptide_info_list):
            self._add_elution_profile(
                fig, row=i+1,col=1,
                tims_data=tims_data, 
                peptide_info=peptide_info, 
                mz_tol=mz_tol, rt_tol=rt_tol,
                im_tol=im_tol,
                include_precursor=include_precursor,
                include_ms1=include_ms1,
            )
        
        fig.update_layout(
            template=self.theme_template,
            # width=width,
            height=self.plot_height,
            hovermode=self.hovermode,
            showlegend=True,
        )
        return fig

