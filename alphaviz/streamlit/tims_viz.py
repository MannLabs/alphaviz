import pandas as pd
import numpy as np
import os
import streamlit as st
import alphabase.psm_reader
from alphatims.bruker import TimsTOF

from alpharaw.viz.df_utils import (
    make_psm_plot_df_for_peptide,
    make_query_plot_df_for_peptide
)

from alpharaw.viz.xic_plot_tims import (
    XIC_Plot_Tims
)

from alpharaw.utils.centroiding import naive_centroid

from alphaviz.streamlit.raw_viz import (
    model_mgr, load_psm, select_psm, psm_plotter
)

xic_plotter = XIC_Plot_Tims()

def show():
    st.write("# Visualization for TimsTOF data")
    tims_file = st.text_input(
        label="TimsTOF file path (*.d or *.hdf in alphatims format)", 
    )

    psm_file_type = st.selectbox(
        label="PSM file type", 
        index=0,
        options=list(alphabase.psm_reader.psm_reader_provider.reader_dict.keys())
    )
    psm_file = st.text_input(
        label="PSM file path",
    )

    if not tims_file or not psm_file or not psm_file_type: return

    tims_data = load_bruker(tims_file)

    raw_name = os.path.split(tims_file)[1]
    if raw_name.endswith(".d"):
        raw_name = raw_name[:-2]
    elif raw_name.endswith(".d.hdf"):
        raw_name = raw_name[:-6]
    elif raw_name.endswith(".hdf"):
        raw_name = raw_name[:-4]
        
    psm_df = load_psm(psm_file, psm_file_type, raw_name)
    psm_df["im"] = psm_df.mobility

    plot_psm = select_psm(psm_df)

    st.write("# Selected PSM to plot")
    st.dataframe(plot_psm)

    use_peptdeep = st.checkbox("Plot predicted")
    if use_peptdeep:
        pred_inten_df = model_mgr.predict_ms2(plot_psm)
    else:
        pred_inten_df = None
    
    ppm = st.number_input(
        "Fragment tolerance (-ppm to +ppm): ", 
        min_value=0., max_value=150., step=0.1, value=20.
    )

    st.markdown("# Plots")

    if st.checkbox("Plot MS2") and len(plot_psm) > 0:

        plot_unmatched = st.checkbox("Plot unmatched peaks", value=False)
        im = plot_psm.mobility.values[0]
        rt_sec = plot_psm.rt.values[0]*60
        mz = plot_psm.precursor_mz.values[0]
        peak_df = tims_data[
            rt_sec-0.1:rt_sec+0.1, 
            im-0.02:im+0.02,
            mz-0.02:mz+0.02, # precursor_mz
            :, # query_mz
        ].sort_values("mz_values")

        peak_mzs, peak_intens = naive_centroid(
            peak_df.mz_values.values,
            peak_df.intensity_values.values,
            centroiding_ppm=ppm
        )

        st.write(f"Number of peaks before centroiding = {len(peak_df)}.")
        st.write(f"Number of peaks after centroiding = {len(peak_mzs)}.")

        plot_df = make_psm_plot_df_for_peptide(
            peak_mzs, peak_intens,
            sequence=plot_psm.sequence.values[0],
            mods=plot_psm.mods.values[0],
            mod_sites=plot_psm.mod_sites.values[0],
            charge=plot_psm.charge.values[0],
            fragment_intensity_df=pred_inten_df,
            ppm=ppm
        )

        fig = psm_plotter.plot(
            plot_df, 
            plot_psm.sequence.values[0], 
            plot_df.modified_sequence.values[0], 
            plot_unmatched_peaks=plot_unmatched
        )

        st.plotly_chart(fig)

    if st.checkbox("Plot XIC") and len(plot_psm) > 0:
        xic_plotter.ppm = ppm

        plot_df = make_query_plot_df_for_peptide(
            sequence=plot_psm.sequence.values[0],
            mods=plot_psm.mods.values[0],
            mod_sites=plot_psm.mod_sites.values[0],
            charge=plot_psm.charge.values[0],
            rt_sec=plot_psm.rt.values[0]*60,
            mobility=plot_psm.mobility.values[0],
            ms_level=2,
            include_precursor_isotopes=False,
            fragment_intensity_df=pred_inten_df,
        )

        st.markdown("### RT view")
        fig = xic_plotter.plot(
            tims_data, 
            query_df=plot_df,
            view_dim="rt",
            title=plot_df.modified_sequence.values[0], 
            add_peak_area=True,
        )

        st.plotly_chart(fig)

        st.markdown("### Mobility view")
        fig = xic_plotter.plot(
            tims_data, 
            query_df=plot_df,
            view_dim="im",
            title=plot_df.modified_sequence.values[0], 
            # add_peak_area=True,
        )

        st.plotly_chart(fig)

        # if st.checkbox("Plot DIA MS1 XIC") and not is_dda:
        #     plot_df = make_query_plot_df_for_peptide(
        #         sequence=plot_psm.sequence.values[0],
        #         mods=plot_psm.mods.values[0],
        #         mod_sites=plot_psm.mod_sites.values[0],
        #         charge=plot_psm.charge.values[0],
        #         rt_sec=plot_psm.rt.values[0]*60,
        #         ms_level=1,
        #         include_precursor_isotopes=True,
        #     )

        #     fig = xic_plotter.plot(
        #         spectrum_df, peak_df, 
        #         query_df=plot_df,
        #         title=plot_df.modified_sequence.values[0], 
        #         add_peak_area=True,
        #     )

        #     st.plotly_chart(fig)

@st.cache_data
def load_bruker(bruker_path):
    return TimsTOF(bruker_path, slice_as_dataframe=True)