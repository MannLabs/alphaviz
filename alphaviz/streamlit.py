import pandas as pd
import streamlit as st
import alphabase.psm_reader
import alpharaw.utils.ms_path_utils
from alpharaw.viz.psm_plot import PSM_Plot
from alpharaw.viz.xic_plot import XIC_Plot
from alpharaw.viz.df_utils import (
    make_psm_plot_df_for_peptide,
    make_xic_plot_df_for_peptide
)
from alpharaw.wrappers.alphatims_wrapper import AlphaTimsWrapper
from peptdeep.pretrained_models import ModelManager

import alphaviz.alphax_utils as xutils

psm_plotter = PSM_Plot()
xic_plotter = XIC_Plot()
model_mgr = ModelManager(mask_modloss=False)

def run():
    raw_file = st.text_input(
        label="Raw file path", 
        value="/Users/wenfengzeng/data/ap_msdata/HeLa_500ng.raw.hdf"
    )
    dda = st.checkbox(label="DDA", value=True)

    psm_file = st.text_input(
        label="PSM file path",
        value="/Users/wenfengzeng/data/ap_msdata/msms.txt"
    )
    psm_file_type = st.selectbox(
        label="PSM file type", 
        index=1,
        options=list(alphabase.psm_reader.psm_reader_provider.reader_dict.keys())
    )
    use_peptdeep = st.checkbox("Plot predicted")

    loaded = st.checkbox("Load data", value=True)

    if loaded:
        raw_name = alpharaw.utils.ms_path_utils.get_raw_name(raw_file)
        msdata = xutils.get_msdata(raw_file)
        tims_data = AlphaTimsWrapper(msdata, dda=dda)

        psm_df = alphabase.psm_reader.psm_reader_provider.get_reader(
            psm_file_type
        ).import_file(psm_file)
        psm_df = psm_df.query(f"raw_name=='{raw_name}'").reset_index(drop=True)

        key_cols = ["sequence","mods","mod_sites","charge","rt"]
        psm_df = pd.concat(
            [
                psm_df[key_cols], 
                psm_df[[col for col in psm_df.columns if col not in key_cols]]
            ],
            axis=1
        )
    else:
        psm_df = pd.DataFrame()
    
    st.dataframe(psm_df)

        # st.data_editor(
        #     psm_df,
        #     hide_index=False,
        #     disabled=psm_df.columns,
        #     num_rows="dynamic",
        # )

    psm_id = st.selectbox(
        "PSM ID to plot", options=psm_df.index.values,
        index=None
    )
    if psm_id is None: return
    psm_id = int(psm_id)
        
    plot_psm = psm_df.loc[psm_id:psm_id]

    st.write("# Selected PSM to plot")
    st.dataframe(plot_psm)

    if st.checkbox("Plot MS2", value=False) and len(plot_psm) > 0:
        if use_peptdeep:
            pred_inten_df = model_mgr.predict_ms2(plot_psm)
        else:
            pred_inten_df = None

        spec_idx = plot_psm.spec_idx.values[0]
        peak_mzs, peak_intens = msdata.get_peaks(spec_idx)

        plot_df = make_psm_plot_df_for_peptide(
            peak_mzs, peak_intens,
            sequence=plot_psm.sequence.values[0],
            mods=plot_psm.mods.values[0],
            mod_sites=plot_psm.mod_sites.values[0],
            charge=plot_psm.charge.values[0],
            fragment_intensity_df=pred_inten_df,
        )

        fig = psm_plotter.plot(
            plot_df, 
            plot_psm.sequence.values[0], 
            plot_df.modified_sequence.values[0], 
            plot_unmatched_peaks=True
        )

        st.plotly_chart(fig)

    if st.checkbox("Plot XIC", value=True) and len(plot_psm) > 0:
        plot_df = make_xic_plot_df_for_peptide(
            sequence=plot_psm.sequence.values[0],
            mods=plot_psm.mods.values[0],
            mod_sites=plot_psm.mod_sites.values[0],
            charge=plot_psm.charge.values[0],
            rt_sec=plot_psm.rt.values[0]*60,
            ms_level=1 if dda else 2,
            include_precursor_isotopes=True if dda else False,
        )

        fig = xic_plotter.plot(
            tims_data=tims_data, 
            query_df=plot_df,
            title=plot_df.modified_sequence.values[0], 
        )

        st.plotly_chart(fig)

run()