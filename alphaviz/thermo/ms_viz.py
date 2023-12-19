import pandas as pd
import numpy as np

import plotly.graph_objects as go

from alphabase.peptide.fragment import (
    get_charged_frag_types
)

from peptdeep.pretrained_models import ModelManager

from alpharaw.match.psm_match_alphatims import PepSpecMatch_AlphaTims

from .ms2_plot import MS2_Plot
from .xic_plot import XIC_1D_Plot
from .reader_utils import load_psms

from .peptdeep_utils import (
    match_ms2, get_frag_df_from_pep_frag_df,
    predict_one_peptide, get_pep_frag_df_from_dfs,
)

class MS_Viz:
    _min_frag_mz:float = 200.0
    _labeled_sites = ['K','N-term']
    remove_unlabeled_fragments = False
    ms1_use_ppm = True
    ms2_use_ppm = True
    ms2_tol=20.0
    ms1_tol=20.0
    rt_sec_tol_to_slice_spectrum = 3.0
    im_tol_to_slice_spectrum = 0.05
    find_closest_ms2_by_rt_sec = True
    use_predicted_ms2 = True
    def __init__(self, 
        model_mgr:ModelManager,
        frag_types:list = ['b','y','b-modloss','y-modloss'],
    ):
        self.model_mgr = model_mgr
        self.tims_data = None
        self.pep_spec_match = None

        self.psm_df = pd.DataFrame()
        self.fragment_mz_df = pd.DataFrame()
        self.fragment_intensity_df = pd.DataFrame()

        self._frag_types = frag_types
        self._max_frag_charge = 2 # fixed
        self.charged_frag_types = get_charged_frag_types(
            self._frag_types, self._max_frag_charge
        )

        self.prediction_mode=False

        self.ms2_plot = MS2_Plot()
        self.xic_1d_plot = XIC_1D_Plot()

    @property
    def min_frag_mz(self):
        return self._min_frag_mz
    
    @min_frag_mz.setter
    def min_frag_mz(self, val):
        self._min_frag_mz = val
        self.xic_1d_plot.min_frag_mz = val

    def load_ms_data(self, ms_file, ms_file_type, dda:bool):
        self.tims_match = PepSpecMatch_AlphaTims(
            self.charged_frag_types, match_closest=True,
        )
        self.tims_match.load_ms_data(ms_file, ms_file_type, dda)
        self.tims_data = self.tims_match.tims_data

    def add_rt_im_columns_to_psm_df_if_missing(self):
        self.psm_df = self.tims_match._add_missing_columns_to_psm_df(
            self.psm_df
        )

    def load_psms(self, 
        psm_file:str, psm_type:str,
        get_fragments:bool=False,
        add_modification_mapping:dict=None,
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
            add_modification_mapping=add_modification_mapping,
        )

    def predict_one_pep_frag_df(self,
        one_pept_df:pd.DataFrame
    )->pd.DataFrame:
        return predict_one_peptide(
            self.model_mgr, one_pept_df, 
            self.tims_data.rt_max_value,
            self.prediction_mode,
            self._labeled_sites if self.remove_unlabeled_fragments else None
        )

    def pep_frag_df_from_dfs(self,
        one_pept_df:pd.DataFrame,
        frag_mz_df: pd.DataFrame,
        frag_inten_df:pd.DataFrame = None,
    ):
        if frag_inten_df is None:
            frag_inten_df = frag_mz_df.copy()
            frag_inten_df.values[:] = 1.0

        return get_pep_frag_df_from_dfs(
            one_pept_df, frag_mz_df,
            frag_inten_df,
            self.tims_data.rt_max_value,
            use_predicted_values=False,
        )

    def extract_one_pep_frag_df(self,
        one_pept_df:pd.DataFrame,
    )->pd.DataFrame:
        return get_pep_frag_df_from_dfs(
            one_pept_df,
            self.fragment_mz_df, 
            self.fragment_intensity_df,
            self.tims_data.rt_max_value,
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
        pep_frag_df: pd.DataFrame,
    ):
        raise NotImplementedError('TODO for timsTOF data')

    def plot_elution_profile(self,
        pep_frag_df: pd.DataFrame,
        include_precursor:bool=True,
        include_ms1:bool=True,
    )->go.Figure:
        """Based on `alphaviz.plotting.plot_elution_profile`

        Parameters
        ----------
        pep_frag_df : pd.DataFrame
            alphaviz pep_frag_df, 
            see `self.predict_one_peptide`.

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
        self.xic_1d_plot.ms1_tol = self.ms1_tol
        self.xic_1d_plot.ms1_use_ppm = self.ms1_use_ppm
        self.xic_1d_plot.ms2_tol = self.ms2_tol
        self.xic_1d_plot.ms2_use_ppm = self.ms2_use_ppm
        return self.xic_1d_plot.plot(
            self.tims_data,
            pep_frag_df=pep_frag_df,
            include_precursor=include_precursor,
            include_ms1=include_ms1,
        )

    def _add_unmatched_df(self, plot_df, spec_df):
        spec_df['ions'] = "-"
        spec_df['fragment_indices'] = -1
        return pd.concat([spec_df, plot_df], ignore_index=True)

    def plot_mirror_ms2(self, 
        pep_frag_df:pd.DataFrame,
        frag_df:pd.DataFrame=None, 
        spec_df:pd.DataFrame=None, 
        title:str="", 
        matching_mode:str="closest",
        plot_unmatched_peaks:bool=False,
    )->go.Figure:
        """Plot mirrored MS2 for PSMs. 
        Top: experimentally matched 

        Parameters
        ----------

        pep_frag_df : pd.DataFrame
            pep_frag_df in alphaviz format

        frag_df : pd.DataFrame, optional
            Fragment DF
        
        spec_df : pd.DataFrame, optional
            AlphaTims sliced DataFrame for raw data,
            by default None

        matching_mode : {'closest', 'highest'}, optional
            peak matching mode, by default "closest"
        
        plot_unmatched_peaks : bool, optional
            by default True

        Returns
        -------
        go.Figure
            plotly Figure object
        """
        if frag_df is None:
            frag_df = get_frag_df_from_pep_frag_df(pep_frag_df)
        if spec_df is None:
            spec_df = self.get_ms2_spec_df(pep_frag_df)

        frag_df = frag_df.query(f'mz_values>={max(spec_df.mz_values.min()-0.1, self._min_frag_mz)}')
        spec_df['intensity_values'] = spec_df.intensity_values.astype(np.float64)
        plot_df, pcc, spc = match_ms2(
            spec_df=spec_df, frag_df=frag_df,
            use_ppm=self.ms2_use_ppm,
            mz_tol=self.ms2_tol, 
            matching_mode=matching_mode,
        )

        pep_frag_df['pcc'] = pcc
        pep_frag_df['spc'] = spc

        if plot_unmatched_peaks:
            plot_df = self._add_unmatched_df(
                plot_df, spec_df
            )

        if not title:
            title = f"{pep_frag_df['mod_seq_charge'].values[0]} PCC={pcc:.3f}"

        plot_df = plot_df.query('intensity_values!=0')

        return self.ms2_plot.plot(
            plot_df, 
            title=title,
            sequence=pep_frag_df['sequence'].values[0],
            plot_unmatched_peaks=plot_unmatched_peaks,
        )

    def get_ms2_spec_df(self, pep_frag_df)->pd.DataFrame:
        im_slice = (
            slice(None) if pep_frag_df['im'].values[0] == 0 else 
            slice(
                pep_frag_df['im'].values[0]-self.im_tol_to_slice_spectrum,
                pep_frag_df['im'].values[0]+self.im_tol_to_slice_spectrum
            )
        )
        query_rt = pep_frag_df['rt_sec'].values[0]
        rt_slice = slice(
            query_rt-self.rt_sec_tol_to_slice_spectrum,
            query_rt+self.rt_sec_tol_to_slice_spectrum
        )

        spec_df = self.tims_data[
            rt_slice, im_slice, 
            pep_frag_df['precursor_mz'].values[0]:pep_frag_df['precursor_mz'].values[0],
        ]
        # spec_df = spec_df[
        #     (spec_df.quad_low_mz_values <= pep_frag_df['precursor_mz'].values[0])
        #     &(spec_df.quad_high_mz_values >= pep_frag_df['precursor_mz'].values[0])
        # ].reset_index(drop=True)

        _df = spec_df

        if self.find_closest_ms2_by_rt_sec:
            min_rt_left_dist = 1000000
            min_rt_right_dist = 1000000
            for _, df in spec_df.groupby('frame_indices'):
                if (
                    abs(df.rt_values.values[0]-query_rt) < min_rt_right_dist 
                    and df.rt_values.values[0]>=query_rt
                ):
                    _df_right = df
                    min_rt_right_dist = abs(df.rt_values.values[0]-query_rt)
                if (
                    abs(df.rt_values.values[0]-query_rt) < min_rt_left_dist 
                    and df.rt_values.values[0]<=query_rt
                ):
                    _df_left = df
                    min_rt_left_dist = abs(df.rt_values.values[0]-query_rt)
            if '_df_left' in locals() and '_df_right' in locals():
                if _df_left is _df_right:
                    _df = _df_left
                else:
                    _df = pd.concat([_df_left, _df_right])
            elif '_df_left' in locals():
                _df = _df_left
            elif '_df_right' in locals():
                _df = _df_right
        return _df

