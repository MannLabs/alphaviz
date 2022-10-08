import torch

import pandas as pd
import numpy as np

from peptdeep.pretrained_models import ModelManager

from alphabase.constants.modification import MOD_MASS
from alphabase.peptide.precursor import calc_precursor_mz

from typing import Union, Tuple

from peptdeep.mass_spec.match import (
    match_centroid_mz, match_profile_mz
)

from peptdeep.model.ms2 import (
    pearson_correlation, spearman_correlation
)

from alphatims.bruker import TimsTOF

def _get_abc_ion_mask(sequence, labeled_sites):
    mask = np.zeros(len(sequence)-1)
    for i,aa in enumerate(sequence[:-1]):
        if aa in labeled_sites: break
        mask[i] = 0
    return mask

def _get_xyz_ion_mask(sequence, labeled_sites):
    mask = np.zeros(len(sequence)-1)
    for i,aa in enumerate(sequence[:0:-1]):
        if aa in labeled_sites: break
        mask[len(sequence)-i-2] = 0
    return mask

def get_unlabeled_fragment_mask(
    sequence:str,
    fragment_df:pd.DataFrame, 
    labeled_sites:list,
):
    masks = np.ones_like(fragment_df.values)
    if 'N-term' in labeled_sites: pass
    else:
        abc_mask = _get_abc_ion_mask(sequence, labeled_sites)
        for i,ion_type in enumerate(fragment_df.columns):
            if ion_type[0] in 'abc':
                masks[:,i] = abc_mask
    if 'C-term' in labeled_sites: pass
    else:
        xyz_mask = _get_xyz_ion_mask(sequence, labeled_sites)
        for i,ion_type in enumerate(fragment_df.columns):
            if ion_type[0] in 'xyz':
                masks[:,i] = xyz_mask
    return masks


def get_peptide_info_from_dfs(
    one_pept_df: Union[pd.DataFrame, pd.Series], 
    fragment_mz_df: pd.DataFrame, 
    fragment_intensity_df:pd.DataFrame, 
    max_rt_in_seconds:float,
    use_predicted_values:bool=False,
)->pd.DataFrame:
    """

    Returns
    -------
    pd.DataFrame
        peptide_info df in alphaviz format
    """
    df = parse_one_pept_df(one_pept_df)
    df = calc_precursor_mz(df)
    frag_mz_df = fragment_mz_df.iloc[
        df.frag_start_idx.values[0]:df.frag_end_idx.values[0],:
    ].reset_index(drop=True)
    frag_inten_df = fragment_intensity_df.iloc[
        df.frag_start_idx.values[0]:df.frag_end_idx.values[0],:
    ].reset_index(drop=True)

    df['mod_seq'] = [get_mod_seq(
        df['sequence'].values[0], 
        df['mods'].values[0],
        df['mod_sites'].values[0],
    )]
    df['mod_seq_charge'] = df['mod_seq'].str.cat(','+df['charge'].astype('U'))

    df["rt_pred"] *= max_rt_in_seconds
    
    if use_predicted_values:
        df['rt_sec'] = df['rt_pred']
        df['im'] = df['mobility_pred']
    else:
        df['rt_sec'] = df.rt*60
        if 'mobility' in df.columns:
            df['im'] = df.mobility
        else:
            df['im'] = 0


    nAA = df["nAA"].values[0]
    charged_frag_types = []
    frag_types = []
    frag_charges = []
    frags = []
    intens = []
    frag_nums = []
    frag_indices = []
    for column in frag_mz_df.columns.values:
        for i,(mz,inten) in enumerate(zip(
            frag_mz_df[column].values,frag_inten_df[column].values
        )):
            if mz < 10: continue
            frag_num = _get_frag_num(column, i, nAA)
            frag_type, charge = column.split('_z')
            charge = int(charge)
            charged_frag_type = (
                frag_type[0] + str(frag_num) +
                frag_type[1:] + '+'*charge
            )
            charged_frag_types.append(charged_frag_type)
            frag_types.append(frag_type)
            frag_charges.append(charge)
            frags.append(mz)
            intens.append(inten)
            frag_nums.append(frag_num)
            frag_indices.append(i)
    df = pd.concat([df]*len(frag_types),ignore_index=True)
    df['frag_type'] = frag_types
    df['frag_charge'] = frag_charges
    df['ion'] = charged_frag_types
    df['frag_mz'] = frags
    df['frag_intensity'] = intens
    df['frag_number'] = frag_nums
    df['frag_index'] = frag_indices
    return df

def predict_one_peptide(
    model_mgr:ModelManager, 
    one_pept_df:Union[pd.DataFrame, pd.Series],
    max_rt_in_seconds,
    use_predicted_values:bool=False,
    labeled_sites:list = None,
)->pd.DataFrame:
    """Predict RT/Mobility/MS2 for one peptide (df)

    Parameters
    ----------
    one_pept_df : pd.DataFrame
        AlphaBase DF containing only one peptide (one row)

    max_rt_in_seconds : float
        use to convert 'rt_pred' predicted from peptdeep 
        to a real RT values (seconds).

    Returns
    -------
    pd.DataFrame
        peptide_info in alphaviz format
    """
    df = parse_one_pept_df(one_pept_df)
    predict_dict = model_mgr.predict_all(
        df, predict_items=['mobility','rt','ms2'],
        multiprocessing=False
    )

    precursor_df = predict_dict["precursor_df"]
    frag_mz_df = predict_dict['fragment_mz_df']
    frag_inten_df = predict_dict['fragment_intensity_df']

    if labeled_sites:
        masks = get_unlabeled_fragment_mask(
            precursor_df.sequence.values[0],
            frag_mz_df, labeled_sites
        )

        frag_mz_df.values[:] = frag_mz_df.values*masks
        frag_inten_df.values[:] = frag_inten_df.values*masks

    return get_peptide_info_from_dfs(
        precursor_df, 
        frag_mz_df, 
        frag_inten_df, 
        max_rt_in_seconds, 
        use_predicted_values
    )

def parse_one_pept_df(
    one_pept_df:Union[pd.Series,pd.DataFrame]
)->pd.DataFrame:
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
    return df 

def _get_frag_num(frag_type, idx, nAA):
    if frag_type[0] in "abc":
        frag_num = idx+1
    else:
        frag_num = nAA-idx-1
    return frag_num

def get_mod_seq(sequence, mods, mod_sites, mod_as_mass=True, **kwargs):
    seq = '_'+sequence+'_'
    if len(mods) == 0: return seq
    mods = mods.split(';')
    if mod_as_mass: 
        mods = [round(MOD_MASS[mod]) for mod in mods]
        mods = [
            '+'+str(mod) if mod > 0 else str(mod) 
            for mod in mods
        ]
    mod_sites = [int(i) for i in mod_sites.split(';')]
    mod_sites = [i if i>=0 else len(seq)-1 for i in mod_sites]
    mods = list(zip(mod_sites, mods))
    mods.sort(reverse=True)

    for i, mod in mods:
        seq = seq[:i+1]+'('+mod+')'+seq[i+1:]
    return seq

def get_frag_df_from_peptide_info(
    peptide_info:pd.DataFrame
)->pd.DataFrame:
    return peptide_info[
        ['frag_intensity', 'frag_mz', 
        'frag_number', 'frag_index',
        'ion']
    ].rename(columns={
        'frag_intensity':'intensity_values',
        'frag_mz':'mz_values',
        'frag_number':'fragment_numbers',
        'frag_index':'fragment_indices',
        'ion':'ions',
    })

def match_tims_ms2_for_tuning(
    ms_data: TimsTOF,
    psm_df: pd.DataFrame,
    frag_mz_df: pd.DataFrame, 
    mz_tol: float=20.0, 
):
    frag_inten_df = pd.DataFrame(
        np.zeros_like(frag_mz_df.values),
        columns=frag_mz_df.columns
    )



def match_ms2(
    spec_df: pd.DataFrame, 
    frag_df: pd.DataFrame, 
    mz_tol=50, 
    matching_mode="centroid",
)->Tuple[pd.DataFrame, float, float]:
    """

    Parameters
    ----------
    spec_df : pd.DataFrame
        AlphaTims DF object

    frag_df : pd.DataFrame
        Fragment DF of peptides, similar to spec_df

    Returns
    -------
    Tuple[pd.DataFrame, float, float]
        pd.DataFrame: DF contains matched, predicted and unmatched peaks
        float: Pearson correlation
        float: Spearman correlation
    """
    frag_df = frag_df.copy()
    spec_df.sort_values('mz_values', inplace=True)
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

    matched_df = frag_df.copy()
    matched_df['mz_values'] = matched_mzs
    matched_df['intensity_values'] = matched_intens
    matched_df['mass_dev_ppm'] = mass_errs
    matched_df = matched_df[matched_idxes!=-1]

    if len(matched_df) == 0:
        return matched_df, 0, 0

    pcc = pearson_correlation(
        torch.tensor(matched_intens.reshape(1,-1)),
        torch.tensor(frag_df.intensity_values.values.reshape(1,-1))
    ).item()

    spc = spearman_correlation(
        torch.tensor(matched_intens.reshape(1,-1)),
        torch.tensor(frag_df.intensity_values.values.reshape(1,-1)),
        device=torch.device('cpu')
    ).item()

    frag_df['intensity_matched'] = matched_intens
    frag_df['intensity_values'] *= (
        -matched_intens.max()/frag_df.intensity_values.max()
    )

    df_list = [matched_df,frag_df]
    plot_df = pd.concat(
        df_list
    ).reset_index(drop=True)

    return plot_df, pcc, spc
