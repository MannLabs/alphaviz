#!python
"""
This module provides pytest tests for the functions from preprocessing.py file
"""

import pytest
import alphaviz.preprocessing


def test_extract_analysis_unique_proteins():
    mq_analysis_file = "../test_data/evidence.txt"
    proteins = alphaviz.preprocessing.extract_analysis_unique_proteins(mq_analysis_file)

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


def test_extract_identified_ions():
    peptide = 'NTINHN'
    all_ions = ['', 'b2-H2O', '', '', 'b2', 'b3', '', 'b5-NH3', '', '', '', '', 'y1-NH3', 'y1', 'y4-H2O', 'y4']
    assert [False,True,True,False,True,False] == alphaviz.preprocessing.extract_identified_ions(all_ions, peptide, 'b'), \
    "The extraction of identified b-ions doesn't work"
    assert [False,False,True,False,False,True] == alphaviz.preprocessing.extract_identified_ions(all_ions, peptide, 'y'), \
    "The extraction of identified y-ions doesn't work"
