import alphabase.psm_reader

import alpharaw.thermo
import alpharaw.sciex
from alpharaw.ms_data_base import ms_reader_provider

import alpharaw.wrappers.alphatims_wrapper

_raw_exts={
    ".raw": "thermo",  
    ".wiff": "sciex", 
    ".raw.hdf": "hdf", 
    ".wiff.hdf": "hdf"
}

def get_msdata(raw_path):
    for _ext, _type in _raw_exts.items():
        if raw_path.endswith(_ext): 
            return ms_reader_provider.get_reader(_type)
    return None

def wrap_raw_to_timstof(raw_path, dda=False):
    msdata = get_msdata(raw_path)
    if msdata is None: return None
    msdata.import_raw(raw_path)
    return alpharaw.wrappers.alphatims_wrapper.AlphaTimsWrapper(msdata, dda=dda)


