##### TODO

import pandas as pd
import numpy as np
import numba
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px

from alphatims.bruker import TimsTOF

from alpharaw.viz.xic_plot_tims import get_plotting_slices

def plot_tims_profile_heatmap(
    tims_data: TimsTOF,
    query_df: dict,
    mz_ppm: float = 20.0,
    rt_sec_win: float = 30.0,
    im_win: float = 0.05,
    n_cols: int = 5,
    # width: int = 180,
):
    img_datas = []
    img_titles = []
    mz = query_df.precursor_mz.values[0]
    rt_slice, im_slice, prec_mz_slice, _ = get_plotting_slices(
        tims_data, 
        rt_sec=query_df.rt_sec.values[0],
        im=query_df.im.values[0],
        rt_sec_win=rt_sec_win,
        im_win=im_win,
        precursor_left_mz=mz*(1-mz_ppm*1e-6),
        precursor_right_mz=mz*(1+mz_ppm*1e-6),
    )

    min_frame, max_frame = np.searchsorted(
        tims_data.rt_values, [rt_slice.start, rt_slice.stop]
    )

    max_scan, min_scan = (
        tims_data.scan_max_index - np.searchsorted(
            tims_data.mobility_values[::-1],
            [im_slice.start, im_slice.stop],
            side="right" 
        )
    )

    # create an elution profile for the precursor
    # tims_df = tims_data[
    #     rt_slice,
    #     im_slice,
    #     0,
    #     prec_mz_slice,
    # ]

    # img_datas.append(_make_dense(tims_df))
    # img_titles.append("precursor")

    # create elution profiles for all fragments
    for ion_name, mz in query_df[["ion_name","mz"]].values:
        mz = float(mz)
        tims_df = tims_data[
            rt_slice,
            im_slice,
            prec_mz_slice,
            slice(mz*(1-mz_ppm*1e-6), mz*(1+mz_ppm*1e-6)),
        ]
        if len(tims_df) > 0:
            img_datas.append(
                _make_dense(
                    tims_df,      
                    min_frame, max_frame,
                    min_scan, max_scan,
                )
            )
            img_titles.append(ion_name)
    
    return px.imshow(
        np.stack(img_datas, axis=0), 
        facet_col=0,facet_col_wrap=n_cols
    )
    
def _make_dense(
    tims_df,
    min_frame, max_frame,
    min_scan, max_scan,
):
    return _make_dense_numba(
        tims_df.frame_indices.values,
        tims_df.scan_indices.values,
        tims_df.intensity_values.values,
        min_frame, max_frame,
        min_scan, max_scan,
    )

@numba.njit
def _make_dense_numba(
    frame_indices,
    scan_indices,
    intensity_values,
    min_frame, max_frame,
    min_scan, max_scan,
):
    ret_arrays = np.zeros(
        (max_frame-min_frame+1,max_scan-min_scan+1),
        dtype=intensity_values.dtype
    )

    for frame, scan, inten in zip(
        frame_indices, scan_indices, intensity_values
    ):
        ret_arrays[frame-min_frame,scan-min_scan] += inten
    
    return ret_arrays