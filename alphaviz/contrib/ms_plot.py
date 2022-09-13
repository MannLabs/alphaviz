import os
from typing import Union

import pandas as pd
import numpy as np
import plotly

import plotly.graph_objects as go

from alphabase.constants.modification import MOD_DF
from alphabase.psm_reader import psm_reader_provider, PSMReaderBase
from alphabase.peptide.fragment import get_charged_frag_types
from alphabase.peptide.fragment import (
    create_fragment_mz_dataframe
)

from peptdeep.pretrained_models import ModelManager
from peptdeep.psm_frag_reader.library_frag_reader import (
    SpectronautMSMSReader
)
from peptdeep.psm_frag_reader.maxquant_frag_reader import (
    MaxQuantMSMSReader
)
from peptdeep.mass_spec.match import (
    match_centroid_mz, match_profile_mz
)

from alphatims.bruker import TimsTOF
from alpharaw.thermo import ThermoRawData
from alpharaw.sciex import SciexWiffData
from alpharaw.wrappers.alphatims_wrapper import AlphaTimsWrapper

from alphaviz.plotting import plot_elution_profile

def get_mod_seq(sequence, mods, mod_sites, **kwargs):
    seq = '_'+sequence+'_'
    if len(mods) == 0: return seq
    mods = mods.split(';')
    mod_sites = [int(i) for i in mod_sites.split(';')]
    mod_sites = [i if i>=0 else len(seq)-1 for i in mod_sites]
    mods = list(zip(mod_sites, mods))
    mods.sort(reverse=True)

    for i, mod in mods:
        seq = seq[:i+1]+'('+mod+')'+seq[i+1:]
    return seq

def add_alphabase_mods(
    psm_reader:PSMReaderBase
):
    for mod in MOD_DF.mod_name.values:
        if mod in psm_reader.modification_mapping:
            psm_reader.modification_mapping[mod].append(
                f'{mod[-1] if mod[-2]=="@" else ""}({mod})'
            )
        else:
            psm_reader.modification_mapping[mod] = [
                f'{mod[-1] if mod[-2]=="@" else ""}({mod})'
            ]
    psm_reader.set_modification_mapping(
        psm_reader.modification_mapping
    )

class MS_Plotter:
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
        self.colorscale_qualitative="Alphabet"
        self.colorscale_sequential="Viridis"

    def load_ms_data(self, ms_file, dda=False):
        if os.path.isfile(ms_file):
            self.ms_data = None
        if ms_file.lower().endswith('.hdf'):
            self.ms_data = TimsTOF(ms_file)
        elif ms_file.lower().endswith('.d'):
            self.ms_data = TimsTOF(ms_file)
        elif ms_file.lower().endswith('.raw'):
            raw_data = ThermoRawData(centroided=True)
            raw_data.import_raw(ms_file)
            self.ms_data = AlphaTimsWrapper(raw_data, dda=dda)
        elif (
            ms_file.lower().endswith('.wiff') or
            ms_file.lower().endswith('.wiff2')
        ):
            raw_data = SciexWiffData(centroided=False)
            raw_data.import_raw(ms_file)
            self.ms_data = AlphaTimsWrapper(raw_data, dda=dda)
        else:
            self.ms_data = None
        
    def load_psms(self, 
        psm_file:str, 
        psm_type:str, 
        get_fragments:bool = False,
        predict_fragments:bool = True,
    )->dict:
        """Load PSMs including DIA PSMs

        Parameters
        ----------
        psm_file : str
            PSM file path

        psm_type : str
            Could be 'alphapept', 'diann', 'maxquant',
            'swath', 'spectronaut', ...
            Both 'swath' and 'spectronaut' use the spec lib (TSV) format.

        get_fragments : bool, optional
            If True, self.fragment_mz_df and 
            self.fragment_intensity_df will also be loaded or predicted.
            By default False
        
        predict_fragments : bool, optional
            If get_fragments is True, the fragments will be 
            loaded from psm_file or predicted.
            If False and the psm_file is MaxQuant msms.txt 
            or spectral library files, it will load intensities from the psm_file.
            Otherwise predict them.
            By default True
        """

        if not get_fragments:
            reader = psm_reader_provider.get_reader(psm_type)
            add_alphabase_mods(reader)
            reader.import_file(psm_file)
            self.psm_df = reader.psm_df
            return

        elif not predict_fragments:
            if os.path.basename(psm_file).lower()=='msms.txt': 
                reader = MaxQuantMSMSReader(
                    frag_types=self._frag_types,
                    max_frag_charge=self._max_frag_charge
                )
            elif psm_type.lower() in ['swath','spectroanut']:
                reader = SpectronautMSMSReader(
                    frag_types=self._frag_types,
                    max_frag_charge=self._max_frag_charge
                )
            else:
                reader = None
            
            if reader is not None:
                reader.import_file(psm_file)
                psm_df = reader.psm_df
                fragment_intensity_df = reader.fragment_intensity_df
                fragment_mz_df = create_fragment_mz_dataframe(
                    psm_df, self.charged_frag_types, 
                    reference_fragment_df=fragment_intensity_df,
                )
                self.psm_df = psm_df
                self.fragment_mz_df = fragment_mz_df
                self.fragment_intensity_df = fragment_intensity_df
                return
        
        reader = psm_reader_provider.get_reader(psm_type)
        add_alphabase_mods(reader)
        reader.import_file(psm_file)
        psm_df = reader.psm_df
        ret = self.model_mgr.predict_all(
            psm_df
        ) 
        self.psm_df = ret['precursor_df']
        self.fragment_mz_df = ret['fragment_mz_df']
        self.fragment_intensity_df = ret['fragment_intensity_df']

    def transfer_learn(self):
        """
        Transfer learning for RT and CCS models based on self.psm_df, 
        and if applicable, MS2 model based on self.fragment_intensity_df
        """
        self.model_mgr.train_ccs_model(
            self.psm_df
        )
        self.model_mgr.train_rt_model(self.psm_df)
        self.model_mgr.train_ms2_model(
            self.psm_df, self.fragment_intensity_df
        )


    def predict_one_peptide(
        self, one_pept_df:pd.DataFrame
    )->dict:
        """Predict RT/Mobility/MS2 for one peptide (df)

        Parameters
        ----------
        one_pept_df : pd.DataFrame
            AlphaBase DF containing only one peptide (one row)

        Returns
        -------
        dict
            peptide_info in alphaviz format
        """
        if len(one_pept_df) == 0: return dict()
        if isinstance(one_pept_df, pd.Series):
            df = pd.DataFrame(
                one_pept_df.to_dict(), index=[0]
            )
        else:
            if len(one_pept_df) > 1:
                df = one_pept_df.iloc[0:1].copy()
            else:
                df = one_pept_df.copy()
        
        predict_dict = self.model_mgr.predict_all(
            df, predict_items=['mobility','rt','ms2'],
            multiprocessing=False
        )
        df = predict_dict["precursor_df"]
        frag_mz_df = predict_dict['fragment_mz_df']
        frag_inten_df = predict_dict['fragment_intensity_df']

        peptide_info = dict(
            sequence = df.sequence.values[0],
            mods = df.mods.values[0],
            mod_sites = df.mod_sites.values[0],
            charge = df.charge.values[0],
        )

        peptide_info['mod_seq'] = get_mod_seq(**peptide_info)
        peptide_info['mod_seq_charge'] = (
            peptide_info['mod_seq'] + '[' + str(peptide_info['charge']) + '+]'
        )

        peptide_info["mz"] = df.precursor_mz.values[0]
        peptide_info["rt_pred"] = df.rt_pred.values[0]
        peptide_info["mobility_pred"] = df.mobility_pred.values[0]
        if 'rt' in df.columns:
            peptide_info['rt_detected'] = df.rt.values[0]
            peptide_info['rt'] = peptide_info['rt_detected']*60
        else:
            peptide_info['rt'] = (
                peptide_info["rt_pred"]*self.ms_data.rt_max_value
            )
        if 'mobility' in df.columns:
            peptide_info['mobility_detected'] = df.mobility.values[0]
            peptide_info['im'] = peptide_info['mobility_detected']
        else:
            peptide_info['im'] = peptide_info['mobility_pred']

        nAA = len(peptide_info["sequence"])
        def get_frag_type(column, idx, nAA):
            if column[0] in "abc":
                ion_num = idx+1
            else:
                ion_num = nAA-idx-1
            charge = int(column[-1])
            return f"{column[0]}{ion_num}{column[1:-3]}{'+'*charge}", ion_num
        frags = {}
        intens = {}
        frag_nums = {}
        for column in frag_mz_df.columns.values:
            for i,(mz,inten) in enumerate(zip(
                frag_mz_df[column].values,frag_inten_df[column].values
            )):
                if mz < 10: continue
                frag_type, ion_num = get_frag_type(column, i, nAA)
                frags[frag_type] = mz
                intens[frag_type] = inten
                frag_nums[frag_type] = ion_num
        peptide_info['fragment_mzs'] = frags
        peptide_info['fragment_intensities'] = intens
        peptide_info['fragment_numbers'] = frag_nums
        peptide_info['fragments'] = peptide_info['fragment_mzs']
        return peptide_info

    def get_frag_df_from_peptide_info(self, peptide_info)->pd.DataFrame:
        return pd.concat([
            pd.DataFrame().from_dict(
                peptide_info['fragment_intensities'], orient='index', 
                columns=['intensity_values']
            ),
            pd.DataFrame().from_dict(
                peptide_info['fragment_mzs'], orient='index', 
                columns=['mz_values']
            ),
            pd.DataFrame().from_dict(
                peptide_info['fragment_numbers'], orient='index', 
                columns=['fragment_numbers']
            ),
        ], axis=1).reset_index().rename(columns={'index':'ions'})

    def plot_elution_profile(
        self,
        peptide_info: dict,
        mz_tol: float = 50,
        rt_tol: float = 30,
        im_tol: float = 0.05,
        # width: int = 900,
        height: int = 400,
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
        return plot_elution_profile(
            self.ms_data, peptide_info=peptide_info,
            mass_dict={}, calculate_fragment_masses=False,
            colorscale_qualitative=self.colorscale_qualitative,
            colorscale_sequential=self.colorscale_sequential,
            mz_tol=mz_tol, rt_tol=rt_tol, im_tol=im_tol,
            title=peptide_info['mod_seq_charge'], 
            height=height,
        )

    def plot_mirror_ms2(self, 
        spec_df:pd.DataFrame, 
        frag_info:Union[pd.DataFrame, dict], 
        title:str="", 
        mz_tol:float=50,
        peak_line_width:float=1.5, 
        height=520,
        color_dict:dict = {
            '-': 'lightgrey', # '-' means umnatched
            'b': 'blue', 'y': 'red',
        },
        template:str="plotly_white", 
        matching_mode:str="centroid",
        plot_unmatched_peaks:bool=True,
    )->go.Figure:
        """Plot mirrored MS2 for PSMs. 
        Top: experimentally matched 

        Parameters
        ----------
        spec_df : pd.DataFrame
            AlphaTims sliced DataFrame for raw data

        frag_info : Union[pd.DataFrame, dict]
            Could be:
            - Fragment DataFrame, see `self.get_frag_df_from_peptide_info`
            - peptide_info, see `self.predict_one_peptide`

        title : str, optional
            figure title, by default ""

        mz_tol : float, optional
            in ppm, by default 50

        peak_line_width : float, optional
            peak line width, by default 1.5

        height : int, optional
            figure height, by default 520

        color_dict : _type_, optional
            Colors for differet peaks, by default 
            { '-': 'lightgrey', # '-' means umnatched 'b': 'blue', 'y': 'red', }

        template : str, optional
            by default "plotly_white"

        matching_mode : str, optional
            peak matching mode, by default "centroid"
        
        plot_unmatched_peaks : bool, optional
            by default True

        Returns
        -------
        go.Figure
            plotly Figure object
        """
        if isinstance(frag_info, pd.DataFrame):
            frag_df = frag_info
        else:
            frag_df = self.get_frag_df_from_peptide_info(frag_info)
            if not title:
                title = frag_info['mod_seq_charge']

        plot_df = self.match_ms2(
            spec_df=spec_df, frag_df=frag_df,
            mz_tol=mz_tol, 
            matching_mode=matching_mode,
            include_unmatched_peak=plot_unmatched_peaks,
        )

        fig = go.Figure()

        if plot_unmatched_peaks:
            self._add_fig_trace(
                fig, plot_df[plot_df.ions=="-"], 
                color_dict['-'], hovertext=False
            )
        self._add_fig_trace(
            fig, plot_df[plot_df.ions.str.startswith('y')],
            color_dict['y'], hovertext=True
        )
        self._add_fig_trace(
            fig, plot_df[plot_df.ions.str.startswith('b')],
            color_dict['b'], hovertext=True
        )

        self._update_fig_layout(
            fig, plot_df, color_dict, template,
            peak_line_width, title, height
        )

        fig_comb = self._add_mass_err_subplot(fig)
        self._add_mass_err_trace(
            fig_comb, plot_df[
                plot_df.ions.str.startswith('b')
            ].query("intensity_values > 0"),
            color_dict['b'],
            name='b ions'
        )
        self._add_mass_err_trace(
            fig_comb, plot_df[
                plot_df.ions.str.startswith('y')
            ].query("intensity_values > 0"),
            color_dict['y'],
            name='y ions'
        )
        fig_comb.update_yaxes(
            title_text=r"ppm", row=5, col=1, 
            range=[-mz_tol, mz_tol]
        )
        return fig_comb
    
    def match_ms2(self, spec_df, frag_df, 
        mz_tol=50, matching_mode="centroid",
        include_unmatched_peak:bool=True,
    ):
        frag_df = frag_df.copy()
        spec_df = spec_df.copy()
        tols = spec_df.mz_values.values*mz_tol*1e-6
        if matching_mode == 'profile':
            matched_idxes = match_profile_mz(
                spec_df.mz_values.values,
                frag_df.mz_values.values, 
                tols
            )
        else:
            matched_idxes = match_centroid_mz(
                spec_df.mz_values.values,
                frag_df.mz_values.values, 
                tols
            )
        matched_intens = spec_df.intensity_values.values[matched_idxes]
        matched_intens[matched_idxes==-1] = 0
        matched_mzs = spec_df.mz_values.values[matched_idxes]
        mass_errs = (
            matched_mzs-frag_df.mz_values.values
        )/frag_df.mz_values.values*1e6
        mass_errs[matched_idxes==-1] = 0
        matched_mzs[matched_idxes==-1] = 0

        frag_df['intensity_matched'] = matched_intens
        frag_df['intensity_values'] *= (
            -matched_intens.max()/frag_df.intensity_values.max()
        )

        matched_df = frag_df.copy()
        matched_df['mz_values'] = matched_mzs
        matched_df['intensity_values'] = matched_intens
        matched_df['mass_dev_ppm'] = mass_errs
        matched_df = matched_df[matched_idxes!=-1]

        df_list = [matched_df,frag_df,]
        if include_unmatched_peak:
            spec_df['ions'] = "-"
            df_list.append(spec_df)
        plot_df = pd.concat(
            df_list
        ).reset_index(drop=True)

        return plot_df

    def _add_mass_err_subplot(self,
        fig,
    ):
        return plotly.subplots.make_subplots(
            rows=6, cols=3, shared_xaxes=True,
            figure=fig,
            specs=[
                [{"rowspan": 4, "colspan": 3}, None, None],
                [None, None, None],
                [None, None, None],
                [None, None, None],
                [{"colspan": 3}, None, None],
                [{}, {}, {}]
            ],
            vertical_spacing=0.07,
            column_widths=[0.25, 0.5, 0.25]
        )

    def _add_mass_err_trace(self,
        subfig, df, color, name, 
    ):
        subfig.add_trace(
            go.Scatter(
                x=df.mz_values,
                y=df.mass_dev_ppm.round(4),
                hovertext=df.ions,
                hovertemplate='<b>m/z:</b> %{x};<br><b>Mass error(ppm):</b> %{y};<br><b>Ion:</b> %{hovertext}.',
                mode='markers',
                marker=dict(
                    color=color,
                    size=6
                ),
                name=name
            ),
            row=5, col=1
        )

    def _update_fig_layout(self,
        fig, plot_df, color_dict, template,
        spectrum_line_width, title, height, 
    ):
        fig.update_layout(
            template=template,
            shapes=[
                dict(
                    type='line',
                    xref='x',
                    yref='y',
                    x0=plot_df.loc[i, 'mz_values'],
                    y0=0,
                    x1=plot_df.loc[i, 'mz_values'],
                    y1=plot_df.loc[i, 'intensity_values'],
                    line=dict(
                        color=color_dict[plot_df.loc[i, 'ions'][0]],
                        width=spectrum_line_width
                    )
                ) for i in plot_df.index
            ],
            xaxis=dict(
                visible=True,
                title='m/z',
            ),
            yaxis=dict(
                title='intensity',
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
            )
        )

    def _add_fig_trace(self,
        fig, df, color, hovertext=True
    ):
        fig.add_trace(
            go.Scatter(
                x=df.mz_values,
                y=df.intensity_values,
                mode='markers',
                marker=dict(color=color, size=1),
                hovertext= df.ions if hovertext else None,
                hovertemplate='<b>m/z:</b> %{x}<br><b>Intensity:</b> %{y}',
                name='',
                showlegend=False
            )
        )


