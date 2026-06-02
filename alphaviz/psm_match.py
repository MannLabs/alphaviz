"""This module provides the `PepSpecMatch_AlphaTims` class for matching peptide-spectrum matches.

Was moved over from alpharaw.
"""
from typing import Tuple, Union

import numpy as np
import pandas as pd
import tqdm
from alphatims.wrapper.alphatims_wrapper import AlphaTimsWrapper

from peptdeep.match.psm_match import PepSpecMatch
from alphatims.bruker import TimsTOF

from alpharaw.ms_data_base import MSData_Base, ms_reader_provider


alphatims_hdf_types = [
    "alphatims",
    "alphatims_hdf",
    "tims.hdf",
]


def load_ms_data_tims(
    ms_file: Union[str, MSData_Base, TimsTOF],
    ms_file_type: str = "alpharaw_hdf",
    dda: bool = False,
    spectra_sorted_by_rt: bool = True,
) -> Tuple[MSData_Base, TimsTOF]:
    """Load ms data as TimsTOF object

    Parameters
    ----------
    ms_file : str
        ms2 file path

    ms_file_type : str, optional
        ms2 file type, could be
        ["alpharaw_hdf","raw.hdf","thermo","sciex","alphapept_hdf","mgf"].
        Default to 'alphatims_hdf'

    dda : bool, optional
        if it is DDA data, by default False

    spectra_sorted_by_rt : bool, optional
        If spectra are already sorted by RT.
        Defaults to True

    Returns
    -------
    tuple
        MSData_Base: alpharaw MS Data (Reader) object
        TimsTOF: AlphaTims object
    """
    if isinstance(ms_file, TimsTOF):
        return None, ms_file
    elif ms_file_type.lower() in alphatims_hdf_types:
        return None, TimsTOF(ms_file)
    else:
        if isinstance(ms_file, MSData_Base):
            raw_data = ms_file
        else:
            raw_data = ms_reader_provider.get_reader(ms_file_type)
            raw_data.import_raw(ms_file)

            if not spectra_sorted_by_rt:
                # RT may not be sorted in AP HDF for timsTOF after preprocessing
                raw_data._sort_rt()

        tims_data = AlphaTimsWrapper(raw_data, dda=dda)
        return raw_data, tims_data


class PepSpecMatch_AlphaTims(PepSpecMatch):
    """
    Inherited from :class:`peptdeep.match.psm_match.PepSpecMatch`, but
    this can be used for DIA PSM matching by selecting
    MS2 spectra with RT (and IM) values.
    """

    #: RT win to get a MS2 spectrum by slicing
    rt_sec_tol_to_slice_ms2 = 3.0

    #: IM win to get a MS2 spectrum by slicing
    im_tol_to_slice_ms2 = 0.03

    #: find closest MS2 for the given RT when slicing
    find_k_nearest_ms2_by_rt = True
    k_rt_nearest = 1

    # : find closest MS2 for the given RT when slicing
    find_k_nearest_ms2_by_im = False
    k_im_nearest = 11

    def get_peak_df(
        self,
        precursor_mz: float,
        rt_sec: float,
        im: float = 0.0,
    ) -> pd.DataFrame:
        """
        Parameters
        ----------
        precursor_mz : float
            Precursor m/z value
        rt_sec : float
            RT value in seconds
        im : float, optional
            Ion mobility, by default 0.0

        Returns
        -------
        pd.DataFrame
            peak_df in alphatims DF format
        """
        rt_slice = slice(
            rt_sec - self.rt_sec_tol_to_slice_ms2,
            rt_sec + self.rt_sec_tol_to_slice_ms2,
        )

        if im == 0 or self.tims_data.scan_max_index == 1:
            im_slice = slice(None)
        elif self.find_k_nearest_ms2_by_im and self.tims_data.scan_max_index > 1:
            # AlphaTims without AlphaRaw for .d files
            im_slice = self.tims_data.scan_max_index - np.searchsorted(
                self.tims_data.mobility_values[::-1], im
            )
        else:
            im_slice = slice(
                im - self.im_tol_to_slice_ms2, im + self.im_tol_to_slice_ms2
            )

        spec_df = self.tims_data[rt_slice, im_slice, precursor_mz:precursor_mz]

        def find_k_nearest(array, val, k=3):
            nearest = np.argmin(np.abs(array - val))
            if nearest <= k // 2:
                return slice(k)
            elif nearest >= len(array) - k // 2 - 1:
                return slice(-k, None)
            else:
                return slice(nearest - k // 2, nearest + k // 2 + 1)

        if (
            self.find_k_nearest_ms2_by_im
            and im > 0
            and self.tims_data.scan_max_index > 1
        ):
            # RAW from AlphaRaw, mobility===0 in AlphaTims wrapper obj
            scan_idxes = np.sort(spec_df.scan_indices.unique())
            if len(scan_idxes) > 1:  # im from psm
                scan_idxes = scan_idxes[
                    find_k_nearest(
                        self.raw_data.spectrum_df.mobility.values[scan_idxes],
                        im,
                        self.k_im_nearest,
                    )
                ]
                spec_df = spec_df[spec_df.scan_indices.isin(scan_idxes)]

        if self.find_k_nearest_ms2_by_rt:
            rt_values = np.sort(spec_df.rt_values.unique())
            if len(rt_values) > 1:
                closest_rts = rt_values[
                    find_k_nearest(rt_values, rt_sec, self.k_rt_nearest)
                ]
                spec_df = spec_df[spec_df.rt_values.isin(closest_rts)]

        return spec_df

    def get_peaks(
        self,
        precursor_mz: float,
        rt_sec: float,
        im: float = 0.0,
    ) -> tuple:
        """
        Parameters
        ----------
        precursor_mz : float
            Precursor m/z value
        rt_sec : float
            RT value in seconds
        im : float, optional
            Ion mobility, by default 0.0

        Returns
        -------
        tuple
            np.ndarray: peak m/z values
            np.ndarray: peak intensity values
        """
        spec_df = self.get_peak_df(precursor_mz, rt_sec, im)
        spec_df = spec_df.sort_values("mz_values").reset_index(drop=True)
        return (spec_df.mz_values.values, spec_df.intensity_values.values)

    def load_ms_data(
        self,
        ms_file,
        ms_file_type,
        dda=False,
        spectra_sorted_by_rt=True,
    ):
        self.raw_data, self.tims_data = load_ms_data_tims(
            ms_file, ms_file_type, dda, spectra_sorted_by_rt
        )

    def match_ms2_one_raw(
        self, psm_df_one_raw: pd.DataFrame, verbose: bool = False
    ) -> tuple:
        """
        Matching psm_df_one_raw against
        self.tims_data and self.raw_data
        after `self.load_ms_data()`

        Parameters
        ----------
        psm_df_one_raw : pd.DataFrame
            psm dataframe
            that contains only one raw file

        Returns
        -------
        tuple:
            pd.DataFrame: psm dataframe with fragment index information.

            pd.DataFrame: fragment mz dataframe.

            pd.DataFrame: matched intensity dataframe.

            pd.DataFrame: matched mass error dataframe.
            np.inf if a fragment is not matched.

        """
        self.psm_df = psm_df_one_raw

        psm_df_one_raw = self._add_missing_columns_to_psm_df(psm_df_one_raw)

        (
            fragment_mz_df,
            matched_intensity_df,
            matched_mz_err_df,
        ) = self._prepare_matching_dfs()

        if (
            "mobility" in psm_df_one_raw.columns
            and "mobility" in self.raw_data.spectrum_df.columns
        ):
            query_columns = [
                "frag_start_idx",
                "frag_stop_idx",
                "precursor_mz",
                "rt",
                "mobility",
            ]
        else:
            query_columns = [
                "frag_start_idx",
                "frag_stop_idx",
                "precursor_mz",
                "rt",
            ]

        psm_iters = psm_df_one_raw[query_columns].values
        if verbose:
            psm_iters = tqdm.tqdm(psm_iters)

        for items in psm_iters:
            frag_start_idx = int(items[0])
            frag_stop_idx = int(items[1])

            spec_mzs, spec_intens = self.get_peaks(
                *items[2:],
            )
            self._match_one_psm(
                spec_mzs,
                spec_intens,
                fragment_mz_df,
                matched_intensity_df,
                matched_mz_err_df,
                frag_start_idx,
                frag_stop_idx,
            )
        return (psm_df_one_raw, fragment_mz_df, matched_intensity_df, matched_mz_err_df)

    def match_ms2_multi_raw(
        self,
        psm_df: pd.DataFrame,
        ms_files: Union[dict, list],
        ms_file_type: str = "alphatims",
        dda: bool = False,
    ):
        """
        Matching PSM dataframe against the ms2 files in ms_files
        This method will store matched values as attributes:
        - self.psm_df
        - self.fragment_mz_df
        - self.matched_intensity_df
        - self.matched_mz_err_df

        Parameters
        ----------
        psm_df : pd.DataFrame
            PSM dataframe

        ms_files : dict | list
            if dict: {raw_name: ms2 path}
            if list: [ms2 path1, ms2 path2]

        ms_file_type : str, optional
            One of ["alphatims_hdf","alpharaw_hdf","thermo","sciex","alphapept_hdf","mgf"]
            Defaults to 'alphapept'.

        Returns
        -------
        tuple:
            pd.DataFrame: psm dataframe with fragment index information.

            pd.DataFrame: fragment mz dataframe.

            pd.DataFrame: matched intensity dataframe.

            pd.DataFrame: matched mass error dataframe.
            np.inf if a fragment is not matched.

        """
        raise NotImplementedError(
            "Not necessary for matching multiple raw files using AlphaTims, "
            "loop through `match_ms2_one_raw()`"
        )
