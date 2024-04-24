import pandas as pd
import numpy as np
import streamlit as st
import alphabase.psm_reader
import alpharaw.utils.ms_path_utils
from alpharaw.viz.psm_plot import PSM_Plot
from alpharaw.viz.xic_plot import XIC_Plot
from alpharaw.viz.df_utils import (
    make_psm_plot_df_for_peptide,
    make_query_plot_df_for_peptide,
)
from alpharaw.match.spec_finder import find_spec_idxes

from peptdeep.pretrained_models import ModelManager
from peptdeep.utils import read_peptide_table

import alphaviz.alphax_utils as xutils

psm_plotter = PSM_Plot()
xic_plotter = XIC_Plot()
model_mgr = ModelManager(mask_modloss=False)


def show():
    st.write("# Visualization for RAW data")
    st.write("Support Thermo, Sciex, etc.")
    raw_file = st.text_input(
        label="Raw file path",
    )
    is_dda = st.checkbox(label="DDA", value=False)

    psm_file_type = st.selectbox(
        label="PSM file type",
        index=0,
        options=list(alphabase.psm_reader.psm_reader_provider.reader_dict.keys())
        + ["tsv", "csv", "txt"],
    )
    if psm_file_type in ["tsv", "csv", "txt"]:
        df = pd.DataFrame(
            {
                "raw_name": ["xxx", "xxx", "xxx"],
                "rt": [1.0, 2.0, 3.0],
                "spec_idx": [1, 2, 3],
                "sequence": ["ACDEFGHIK", "LMNPQRSTVK", "WYVSTR"],
                "mods": ["Carbamidomethyl@C", "Acetyl@Protein N-term;Phospho@S", ""],
                "mod_sites": ["2", "0;7", ""],
                "charge": [2, 3, 1],
            }
        )
        infile_expander = st.expander("TSV/CSV/TXT file example:")
        with infile_expander:
            st.write("`raw_name` is the raw file name without .raw.")
            st.write("`rt` is in minutes, required by DIA data.")
            st.write("`spec_idx` is `scan number`-1, required by DDA data.")
            st.dataframe(df, hide_index=True)

    psm_file = st.text_input(
        label="PSM file path",
    )

    if not raw_file or not psm_file or not psm_file_type:
        return

    raw_name = alpharaw.utils.ms_path_utils.get_raw_name(raw_file)
    spectrum_df, peak_df = load_raw(raw_file)

    psm_df = load_psm(psm_file, psm_file_type, raw_name)

    plot_psm = select_psm(psm_df)

    if plot_psm is None or len(plot_psm) == 0:
        return

    if not is_dda or ("spec_idx" not in psm_df.columns and "rt" in psm_df.columns):
        spec_idxes = find_spec_idxes(
            spectrum_df.rt.values,
            spectrum_df.isolation_lower_mz.values,
            spectrum_df.isolation_upper_mz.values,
            plot_psm.rt.values[0] - xic_plotter.rt_sec_win / 120,
            plot_psm.rt.values[0] + xic_plotter.rt_sec_win / 120,
            plot_psm.precursor_mz.values[0],
            plot_psm.precursor_mz.values[0],
        )

        plot_psm["spec_idx"] = spec_idxes[
            np.abs(spectrum_df.rt.values[spec_idxes] - plot_psm.rt.values[0]).argmin()
        ]

    st.write("# Selected PSM to plot")
    st.dataframe(plot_psm)

    use_peptdeep = st.checkbox("Plot predicted")
    if use_peptdeep:
        pred_inten_df = model_mgr.predict_ms2(plot_psm)
    else:
        pred_inten_df = None

    plot_unmatched = st.checkbox("Plot unmatched peaks")

    if st.checkbox("Plot MS2") and len(plot_psm) > 0:
        spec_idx = plot_psm.spec_idx.values[0]
        peak_mzs, peak_intens = get_peaks(spectrum_df, peak_df, spec_idx)

        plot_df = make_psm_plot_df_for_peptide(
            peak_mzs,
            peak_intens,
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
            plot_unmatched_peaks=plot_unmatched,
        )

        st.plotly_chart(fig)

    if st.checkbox("Plot XIC") and len(plot_psm) > 0:

        if st.button("Reset Retention Time Window (60 seconds):"):
            xic_plotter.rt_sec_win = 60

        rt_sec = plot_psm.rt.values[0] * 60
        rt_range = st.slider(
            "Retention Time Range (in seconds)",
            min_value=rt_sec - 200,
            max_value=rt_sec + 200,
            value=(
                rt_sec - xic_plotter.rt_sec_win / 2,
                rt_sec + xic_plotter.rt_sec_win / 2,
            ),
        )
        plot_df = make_query_plot_df_for_peptide(
            sequence=plot_psm.sequence.values[0],
            mods=plot_psm.mods.values[0],
            mod_sites=plot_psm.mod_sites.values[0],
            charge=plot_psm.charge.values[0],
            rt_sec=np.mean(rt_range),
            ms_level=1 if is_dda else 2,
            include_precursor_isotopes=True if is_dda else False,
            fragment_intensity_df=pred_inten_df,
        )
        xic_plotter.rt_sec_win = rt_range[1] - rt_range[0]
        fig = xic_plotter.plot(
            spectrum_df,
            peak_df,
            query_df=plot_df,
            title=plot_df.modified_sequence.values[0],
            add_peak_area=True,
        )

        st.plotly_chart(fig)

        if st.checkbox("Plot DIA MS1 XIC") and not is_dda:
            plot_df = make_query_plot_df_for_peptide(
                sequence=plot_psm.sequence.values[0],
                mods=plot_psm.mods.values[0],
                mod_sites=plot_psm.mod_sites.values[0],
                charge=plot_psm.charge.values[0],
                rt_sec=plot_psm.rt.values[0] * 60,
                ms_level=1,
                include_precursor_isotopes=True,
            )

            fig = xic_plotter.plot(
                spectrum_df,
                peak_df,
                query_df=plot_df,
                title=plot_df.modified_sequence.values[0],
                add_peak_area=True,
            )

            st.plotly_chart(fig)


@st.cache_data(max_entries=1)
def load_raw(raw_file):
    msdata = xutils.get_msdata(raw_file)
    return msdata.spectrum_df, msdata.peak_df


@st.cache_data(max_entries=3)
def load_psm(psm_file, psm_file_type, raw_name):
    if psm_file_type in ["tsv", "csv", "txt"]:
        psm_df = read_peptide_table(psm_file)
    else:
        psm_df = alphabase.psm_reader.psm_reader_provider.get_reader(
            psm_file_type
        ).import_file(psm_file)
    psm_df = psm_df.query(f"raw_name=='{raw_name}'").reset_index(drop=True)

    key_cols = ["sequence", "mods", "mod_sites", "charge", "rt"]
    psm_df = pd.concat(
        [
            psm_df[key_cols],
            psm_df[[col for col in psm_df.columns if col not in key_cols]],
        ],
        axis=1,
    )
    return psm_df


def select_psm(psm_df):
    protein = st.selectbox(
        "Select protein name",
        options=[""] + list(dict.fromkeys(psm_df.proteins.values)),
        index=None,
    )
    if "uniprot_ids" in psm_df.columns:
        protein_id = st.selectbox(
            "Select protein ID or Uniprot ID",
            options=[""] + list(dict.fromkeys(psm_df.uniprot_ids.values)),
            index=None,
        )
    else:
        protein_id = None

    gene = st.selectbox(
        "Select gene name",
        options=[""] + list(dict.fromkeys(psm_df.genes.values)),
        index=None,
    )

    peptide = st.selectbox(
        "Select peptide sequence",
        options=[""] + list(dict.fromkeys(psm_df.sequence.values)),
        index=None,
    )

    if peptide:
        psm_df_to_display = psm_df.query(f"sequence=='{peptide}'")
    else:
        psm_df_to_display = psm_df

    if protein:
        psm_df_to_display = psm_df_to_display[
            psm_df_to_display.proteins.str.contains(protein)
        ]
    if gene:
        psm_df_to_display = psm_df_to_display[
            psm_df_to_display.genes.str.contains(gene)
        ]
    if protein_id:
        psm_df_to_display = psm_df_to_display[
            psm_df_to_display.uniprot_ids.str.contains(protein_id)
        ]

    st.dataframe(psm_df_to_display)

    psm_id = st.selectbox(
        "Select PSM ID to plot", options=psm_df_to_display.index.values, index=None
    )
    if psm_id is None:
        return
    psm_id = int(psm_id)

    plot_psm = psm_df.loc[psm_id:psm_id]

    return plot_psm

def get_peaks(spec_df, peak_df, spec_idx):
    start, stop = spec_df[["peak_start_idx", "peak_stop_idx"]].values[spec_idx]
    return (peak_df.mz.values[start:stop], peak_df.intensity.values[start:stop])
