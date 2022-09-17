import torch

import pandas as pd

from peptdeep.pretrained_models import ModelManager

from typing import Union, Tuple

from peptdeep.mass_spec.match import (
    match_centroid_mz, match_profile_mz
)

from peptdeep.model.ms2 import (
    pearson_correlation, spearman_correlation
)

def get_peptide_info_from_dfs(
    one_pept_df: Union[pd.DataFrame, pd.Series], 
    fragment_mz_df: pd.DataFrame, 
    fragment_intensity_df:pd.DataFrame, 
    max_rt_in_seconds:float,
)->dict:
    """

    Returns
    -------
    dict
        peptide_info in alphaviz format
    """
    df = parse_one_pept_df(one_pept_df)
    frag_mz_df = fragment_mz_df.iloc[
        df.frag_start_idx.values[0]:df.frag_end_idx.values[0],:
    ].reset_index(drop=True)
    frag_inten_df = fragment_intensity_df.iloc[
        df.frag_start_idx.values[0]:df.frag_end_idx.values[0],:
    ].reset_index(drop=True)

    peptide_info = dict(
        sequence = one_pept_df.sequence.values[0],
        mods = one_pept_df.mods.values[0],
        mod_sites = one_pept_df.mod_sites.values[0],
        charge = one_pept_df.charge.values[0],
    )

    peptide_info['mod_seq'] = get_mod_seq(**peptide_info)
    peptide_info['mod_seq_charge'] = (
        peptide_info['mod_seq'] + ',' + 
        str(peptide_info['charge']) + '+'
    )

    peptide_info["mz"] = df.precursor_mz.values[0]
    peptide_info["rt_pred"] = df.rt_pred.values[0]*max_rt_in_seconds
    peptide_info["rt_norm_pred"] = df.rt_norm_pred.values[0]
    peptide_info["mobility_pred"] = df.mobility_pred.values[0]
    if 'rt' in df.columns:
        peptide_info['rt_detected'] = df.rt.values[0]*60
        peptide_info['rt'] = peptide_info['rt_detected']
    else:
        peptide_info['rt'] = peptide_info['rt_pred']
        
    if 'mobility' in df.columns:
        peptide_info['mobility_detected'] = df.mobility.values[0]
        peptide_info['im'] = peptide_info['mobility_detected']
    else:
        peptide_info['im'] = peptide_info['mobility_pred']

    nAA = len(peptide_info["sequence"])
    frags = {}
    intens = {}
    frag_nums = {}
    frag_indices = {}
    for column in frag_mz_df.columns.values:
        for i,(mz,inten) in enumerate(zip(
            frag_mz_df[column].values,frag_inten_df[column].values
        )):
            if mz < 10: continue
            frag_num = _get_frag_num(column, i, nAA)
            frag_type, charge = column.split('_z')
            frag_type = (
                frag_type[0] + str(frag_num) +
                frag_type[1:] + '+'*int(charge)
            )
            frags[frag_type] = mz
            intens[frag_type] = inten
            frag_nums[frag_type] = frag_num
            frag_indices[frag_type] = i
    peptide_info['fragment_mzs'] = frags
    peptide_info['fragment_intensities'] = intens
    peptide_info['fragment_numbers'] = frag_nums
    peptide_info['fragment_indices'] = frag_indices
    peptide_info['fragments'] = peptide_info['fragment_mzs']
    return peptide_info

def predict_one_peptide(
    model_mgr:ModelManager, 
    one_pept_df:Union[pd.DataFrame, pd.Series],
    max_rt_in_seconds
)->dict:
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
    dict
        peptide_info in alphaviz format
    """
    df = parse_one_pept_df(one_pept_df)
    
    predict_dict = model_mgr.predict_all(
        df, predict_items=['mobility','rt','ms2'],
        multiprocessing=False
    )

    return get_peptide_info_from_dfs(
        predict_dict["precursor_df"], 
        predict_dict['fragment_mz_df'], 
        predict_dict['fragment_intensity_df'], 
        max_rt_in_seconds
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

def get_frag_df_from_peptide_info(
    peptide_info
)->pd.DataFrame:
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
        pd.DataFrame().from_dict(
            peptide_info['fragment_indices'], orient='index', 
            columns=['fragment_indices']
        ),
    ], axis=1).reset_index().rename(columns={'index':'ions'})

def match_ms2(
    spec_df: pd.DataFrame, 
    frag_df: pd.DataFrame, 
    mz_tol=50, matching_mode="centroid",
    include_unmatched_peak:bool=True,
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

    matched_df = frag_df.copy()
    matched_df['mz_values'] = matched_mzs
    matched_df['intensity_values'] = matched_intens
    matched_df['mass_dev_ppm'] = mass_errs
    matched_df = matched_df[matched_idxes!=-1]

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
    if include_unmatched_peak:
        spec_df['ions'] = "-"
        df_list.append(spec_df)
    plot_df = pd.concat(
        df_list
    ).reset_index(drop=True)

    return plot_df, pcc, spc
