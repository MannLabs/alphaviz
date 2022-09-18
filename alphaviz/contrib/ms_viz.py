import pandas as pd

import plotly.graph_objects as go

from alphabase.peptide.fragment import (
    get_charged_frag_types
)

from peptdeep.pretrained_models import ModelManager

from .ms2_plot import MS2_Plot
from .xic_plot import XIC_1D_Plot
from .reader_utils import load_ms_data, load_psms

from .peptdeep_utils import (
    match_ms2, get_frag_df_from_peptide_info,
    predict_one_peptide, get_peptide_info_from_dfs,
)

class MS_Viz:
    def __init__(self, 
        model_mgr:ModelManager,
        frag_types:list = ['b','y','b-modloss','y-modloss'],
    ):
        self.model_mgr = model_mgr
        self.ms_data = None
        self.psm_df = pd.DataFrame()
        self.fragment_mz_df = pd.DataFrame()
        self.fragment_intensity_df = pd.DataFrame()

        self._frag_types = frag_types
        self._max_frag_charge = 2 # fixed
        self.charged_frag_types = get_charged_frag_types(
            self._frag_types, self._max_frag_charge
        )

        self.ms2_plot = MS2_Plot()
        self.xic_1d_plot = XIC_1D_Plot()

    def load_ms_data(self, ms_file, dda:bool):
        self.ms_data = load_ms_data(ms_file, dda=dda)

    def load_psms(self, 
        psm_file, psm_type,
        get_fragments=False,
    ):
        (
            self.psm_df, self.fragment_mz_df, 
            self.fragment_intensity_df
        ) = load_psms(
            psm_file, psm_type, 
            get_fragments=get_fragments,
            model_mgr=self.model_mgr,
            frag_types=self._frag_types,
            max_frag_charge=self._max_frag_charge,
        )

    def predict_one_peptide_info(self,
        one_pept_df:pd.DataFrame
    )->dict:
        return predict_one_peptide(
            self.model_mgr, one_pept_df, 
            self.ms_data.rt_max_value
        )

    def extract_one_peptide_info(self,
        one_pept_df:pd.DataFrame,
    )->dict:
        return get_peptide_info_from_dfs(
            one_pept_df,
            self.fragment_mz_df, 
            self.fragment_intensity_df,
            self.ms_data.rt_max_value,
        )

    def transfer_learn(self):
        """
        Transfer learning for RT and CCS models based on self.psm_df, 
        and if applicable, MS2 model based on self.fragment_intensity_df
        """
        self.model_mgr.train_ccs_model(self.psm_df)
        self.model_mgr.train_rt_model(self.psm_df)
        if len(self.fragment_intensity_df) > 0:
            self.model_mgr.train_ms2_model(
                self.psm_df, self.fragment_intensity_df
            )

    def plot_elution_profile_heatmap(self,
        peptide_info: dict,
        mz_tol: float = 50,
        rt_tol: float = 30,
        im_tol: float = 0.05,
    ):
        raise NotImplementedError('TODO for timsTOF data')

    def plot_elution_profile(self,
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
            alphaviz peptide_info dict, 
            see `self.predict_one_peptide`.

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
        return self.xic_1d_plot.plot(
            self.ms_data,
            peptide_info=peptide_info,
            mz_tol=mz_tol,
            rt_tol=rt_tol,
            im_tol=im_tol,
            include_precursor=include_precursor
        )

    def plot_mirror_ms2(self, 
        peptide_info:dict,
        frag_df:pd.DataFrame=None, 
        spec_df:pd.DataFrame=None, 
        title:str="", 
        mz_tol:float=50,
        matching_mode:str="centroid",
        plot_unmatched_peaks:bool=False,
    )->go.Figure:
        """Plot mirrored MS2 for PSMs. 
        Top: experimentally matched 

        Parameters
        ----------

        peptide_info : dict
            peptide_info dict in alphaviz format

        frag_df : pd.DataFrame, optional
            Fragment DF
        
        spec_df : pd.DataFrame, optional
            AlphaTims sliced DataFrame for raw data,
            by default None

        mz_tol : float, optional
            in ppm, by default 50

        matching_mode : str, optional
            peak matching mode, by default "centroid"
        
        plot_unmatched_peaks : bool, optional
            by default True

        Returns
        -------
        go.Figure
            plotly Figure object
        """
        if frag_df is None:
            frag_df = get_frag_df_from_peptide_info(peptide_info)
        if spec_df is None:
            spec_df = self.get_ms2_spec_df(peptide_info)

        plot_df, pcc, spc = match_ms2(
            spec_df=spec_df, frag_df=frag_df,
            mz_tol=mz_tol, 
            matching_mode=matching_mode,
            include_unmatched_peak=plot_unmatched_peaks,
        )

        if not title:
            title = f"{peptide_info['mod_seq_charge']} PCC={pcc:.3f}"

        plot_df = plot_df.query('intensity_values!=0')

        return self.ms2_plot.plot(
            plot_df, 
            title=title,
            sequence=peptide_info['sequence'],
            plot_unmatched_peaks=plot_unmatched_peaks,
        )

    def get_ms2_spec_df(self, peptide_info)->pd.DataFrame:
        im_slice = (
            slice(None) if peptide_info['im'] == 0 else 
            slice(peptide_info['im']-0.05,peptide_info['im']+0.05)
        )
        rt_slice = slice(peptide_info['rt']-0.5,peptide_info['rt']+0.5)

        spec_df = self.ms_data[
            rt_slice, im_slice
        ]
        return spec_df[
            (spec_df.quad_low_mz_values <= peptide_info['mz'])
            &(spec_df.quad_high_mz_values >= peptide_info['mz'])
        ].reset_index(drop=True)

