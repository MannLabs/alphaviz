#!python
"""
This module provides pytest tests for the functions from preprocessing.py file
"""

import pytest
import alphaviz.preprocessing as preproc


def test_get_mq_unique_proteins():
    mq_analysis_file = "../test_data/evidence.txt"
    proteins = preproc.get_mq_unique_proteins(mq_analysis_file)

    assert len(proteins) == 6, \
        "The number of extracted proteins is wrong."
    assert 'Plectin' in proteins, \
        "A unique protein is absent in the extracted list of proteins."


# def test_preprocess_ckg_output():
#     ckg_output_string_correct = "~Q92934;~Q15149"
#     ckg_output_string_correct_2 = " ~Q92934; ~Q15149"
#     ckg_output_string_wrong = "Q92934; Q15149"
#
#     proteins_correct_input = alphaviz.preprocessing.preprocess_ckg_output(ckg_output_string_correct)
#     assert len(proteins_correct_input) == 2, \
#         "The number of extracted proteins is wrong."
#     assert 'Q92934' in proteins_correct_input, \
#         "A unique protein is absent in the extracted list of proteins."
#
#     assert proteins_correct_input == alphaviz.preprocessing.preprocess_ckg_output(ckg_output_string_correct_2), \
#         "Spaces in the ckg string don't influence on the result of the output."
#
#     with pytest.raises(ValueError):
#         alphaviz.preprocessing.preprocess_ckg_output(ckg_output_string_wrong)


def test_get_identified_ions():
    peptide = 'NTINHN'
    all_ions = ['', 'b2-H2O', '', '', 'b2', 'b3', '', 'b5-NH3', '', '', '', '', 'y1-NH3', 'y1', 'y4-H2O', 'y4']
    assert [False,True,True,False,True,False] == preproc.get_identified_ions(all_ions, peptide, 'b'), \
    "The extraction of identified b-ions doesn't work"
    assert [False,False,True,False,False,True] == preproc.get_identified_ions(all_ions, peptide, 'y'), \
    "The extraction of identified y-ions doesn't work"

def test_convert_diann_mq_mod():
    seq1 = 'VSHGSSPSLLEALSSDFLAC(UniMod:4)K'
    assert 'VSHGSSPSLLEALSSDFLAC[Carbamidomethyl (C)]K' == preproc.convert_diann_mq_mod(seq1)
    seq2 = 'VSVINTVDTSHEDMIHDAQM(UniMod:35)DYYGTR'
    assert 'VSVINTVDTSHEDMIHDAQM[Oxidation (M)]DYYGTR' == preproc.convert_diann_mq_mod(seq2)
    seq3 = 'HAEMPVHTGLK(UniMod:2)'
    assert 'HAEMPVHTGLK[Amidated (C-term)]' == preproc.convert_diann_mq_mod(seq3)
    seq4 = 'HAEMPVHTGLKS(UniMod:23)A'
    assert 'HAEMPVHTGLKS[Dehydrated (ST)]A' == preproc.convert_diann_mq_mod(seq4)
    seq5 = 'HAEMPVHTGLKY(UniMod:23)A'
    assert 'HAEMPVHTGLKY[Dehydrated (Y)]A' == preproc.convert_diann_mq_mod(seq5)

    seq1_several_dif_mods = '(UniMod:1)VSHGSSPSLLEALSSDFLAC(UniMod:4)K'
    assert '[Acetyl (N-term)]VSHGSSPSLLEALSSDFLAC[Carbamidomethyl (C)]K' == preproc.convert_diann_mq_mod(seq1_several_dif_mods)
    seq2_several_same_mods = 'CAALVATAEENLC(UniMod:4)C(UniMod:4)EELSSK'
    assert 'CAALVATAEENLC[Carbamidomethyl (C)]C[Carbamidomethyl (C)]EELSSK' == preproc.convert_diann_mq_mod(seq2_several_same_mods)
    seq_no_mod = 'CVNTTLQIK'
    assert "CVNTTLQIK" == preproc.convert_diann_mq_mod(seq_no_mod)


def test_convert_diann_ap_mod():
    seq1 = '(UniMod:27)ELNMIIMLPDETTDLR'
    assert preproc.convert_diann_ap_mod(seq1) == 'pgELNMIIMLPDETTDLR', "The convertion of the N-term modifications doesn't work."

    seq2 = 'ADFSGM(UniMod:35)SQTDLSLSK'
    assert preproc.convert_diann_ap_mod(seq2) == 'ADFSGoxMSQTDLSLSK', "The convertion of the modifications doesn't work."

    seq3 = 'YYYDGDMIC(UniMod:4)'
    assert preproc.convert_diann_ap_mod(seq3) == 'YYYDGDMIcC', \
    "The convertion of the C-term modifications doesn't work."

    seq1_several_dif_mods = '(UniMod:1)VSHGSSPSLLEALSSDFLAC(UniMod:4)K'
    assert 'aVSHGSSPSLLEALSSDFLAcCK' == preproc.convert_diann_ap_mod(seq1_several_dif_mods)
    seq2_several_same_mods = 'CAALVATAEENLC(UniMod:4)C(UniMod:4)EELSSK'
    assert 'CAALVATAEENLcCcCEELSSK' == preproc.convert_diann_ap_mod(seq2_several_same_mods)
    seq_no_mod = 'CVNTTLQIK'
    assert "CVNTTLQIK" == preproc.convert_diann_ap_mod(seq_no_mod)
