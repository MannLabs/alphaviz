import alpharaw.thermo
import alpharaw.sciex
import alpharaw.wrappers.alphapept_wrapper
from alpharaw.ms_data_base import ms_reader_provider

_raw_exts={
    ".raw": "thermo",  
    ".wiff": "sciex", 
    ".raw.hdf": "hdf", 
    ".wiff.hdf": "hdf",
    ".ms_data.hdf": "alphapept_hdf",
}

def get_msdata(raw_path):
    for _ext, _type in _raw_exts.items():
        if raw_path.endswith(_ext): 
            msdata = ms_reader_provider.get_reader(_type)
            msdata.import_raw(raw_path)
            return msdata
    return None


