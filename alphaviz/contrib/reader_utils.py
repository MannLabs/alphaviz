import os

import pandas as pd

from typing import Union, Tuple

from alphabase.constants.modification import MOD_DF
from alphabase.psm_reader import psm_reader_provider, PSMReaderBase
from alphabase.peptide.fragment import (
    get_charged_frag_types,
    create_fragment_mz_dataframe
)

import alpharaw.thermo
import alpharaw.sciex
import alpharaw.wrappers.alphatims_wrapper
import alpharaw.legacy_msdata.mgf

from alpharaw.ms_data_base import ms_reader_provider

from peptdeep.psm_frag_reader.library_frag_reader import (
    SpectronautMSMSReader
)
from peptdeep.psm_frag_reader.maxquant_frag_reader import (
    MaxQuantMSMSReader
)
from peptdeep.pretrained_models import ModelManager

def load_psms(
    psm_file:str, 
    psm_type:str, 
    get_fragments:bool = False,
    model_mgr:ModelManager=None,
    frag_types:list = ['b','y'],
    max_frag_charge:int = 2,
    add_modification_mapping:dict=None,
)->Tuple[
    pd.DataFrame, pd.DataFrame, 
    pd.DataFrame
]:
    """Load PSMs including DIA PSMs.

    If `get_fragments` is True, there are two ways to get fragments:
        1. if `model_mgr` is None
            - if the `psm_file` contains fragment information, e.g.:
              * MaxQuant msms.txt 
              * spectral library files (TSV) 
              fragment intensities will be loaded from the psm_file
            - otherwise return empty fragment DFs
        2. otherwise fragments will be predicted by `model_mgr`

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
    
    model_mgr : ModelManager, optional
        By default None

    Returns
    -------
    tuple
        pd.DataFrame: psm_df
        pd.DataFrame: fragment_mz_df
        pd.DataFrame: fragment_intensity_df
    """
    if not get_fragments:
        reader = psm_reader_provider.get_reader(psm_type)
        _add_alphabase_mods(reader)
        if reader.add_modification_mapping:
            reader.add_modification_mapping(
                add_modification_mapping
            )
        reader.import_file(psm_file)
        return (
            reader.psm_df, pd.DataFrame(), pd.DataFrame()
        )
    else:
        reader = _get_psm_reader(
            psm_file, psm_type, 
            frag_types, max_frag_charge,
            model_mgr is None,
        )
        reader.import_file(psm_file)
        psm_df = reader.psm_df
        if model_mgr is None:
            frag_inten_df = reader.fragment_intensity_df
            frag_mz_df = create_fragment_mz_dataframe(
                psm_df, get_charged_frag_types(
                    frag_types, max_frag_charge
                ),
                reference_fragment_df=frag_inten_df,
            )
            return (
                psm_df, frag_mz_df,
                reader.frag_inten_df
            )
        else:
            ret = model_mgr.predict_all(reader.psm_df) 
            return (
                ret['precursor_df'],
                ret['fragment_mz_df'],
                ret['fragment_intensity_df']
            )

def _add_alphabase_mods(
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

def _get_psm_reader(
    psm_file, psm_type, 
    frag_types, max_frag_charge,
    use_msms_reader
):
    if use_msms_reader:
        if os.path.basename(psm_file).lower()=='msms.txt': 
            reader = MaxQuantMSMSReader(
                frag_types=frag_types,
                max_frag_charge=max_frag_charge
            )
        elif psm_type.lower() in ['swath','spectroanut']:
            reader = SpectronautMSMSReader(
                frag_types=frag_types,
                max_frag_charge=max_frag_charge
            )
    else:
        reader = psm_reader_provider.get_reader(psm_type)
    _add_alphabase_mods(reader)
    return reader