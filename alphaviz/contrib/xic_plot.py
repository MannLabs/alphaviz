import os
import torch
import plotly

import pandas as pd
import numpy as np

import plotly.graph_objects as go

from alphatims.bruker import TimsTOF

from alphaviz.plotting import (
    plot_elution_profile,
    plot_elution_profile_heatmap,
)

class XIC_1D_Plot():
    # hovermode = "x" | "y" | "closest" | False | "x unified" | "y unified"
    hovermode = 'closest'
    plot_height = 550
    colorscale_qualitative="Alphabet"
    colorscale_sequential="Viridis"
    theme_template='plotly_white'

    def plot(self,
        ms_data:TimsTOF,
        peptide_info: dict,
        mz_tol: float = 50,
        rt_tol: float = 30,
        im_tol: float = 0.05,
        include_precursor:bool=True,
    )->go.Figure:
        """Based on `alphaviz.plotting.plot_elution_profile`

        Parameters
        ----------
        peptide_info : dict
            alphaviz peptide_info dict

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
        fig = plot_elution_profile(
            ms_data, peptide_info=peptide_info,
            mass_dict={}, calculate_fragment_masses=False,
            colorscale_qualitative=self.colorscale_qualitative,
            colorscale_sequential=self.colorscale_sequential,
            mz_tol=mz_tol, rt_tol=rt_tol, im_tol=im_tol,
            title=peptide_info['mod_seq_charge'], 
            height=self.plot_height,
            hovermode=self.hovermode,
            theme_template=self.theme_template,
            include_precursor=include_precursor,
        )
        fig.add_vline(
            peptide_info['rt']/60, line_dash="dash", 
            line_color="grey"
        )
        return fig
