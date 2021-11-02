#!python
"""
This module provides pytest tests for the functions from io.py file
"""

import alphaviz.io

# test dataset
mq_evidence_file = "./test_data/evidence.txt" #change this file


def test_read_file():
    mq_columns = ['Protein names', 'Raw file']
    data = alphaviz.io.read_file(mq_analysis_file, mq_columns)

    assert data.shape == (16242, 2), \
        "The number of columns/rows in the read dataset is wrong."
    assert data.columns.tolist() == mq_columns, \
        "The wrong column names were extracted from the dataset."
    assert data['Protein names'].nunique() == 6, \
        "The data of the extracted column are wrong."


def test_upload_evidence_file():
    raw_file_name = 'raw_0'
    data = alphaviz.io.upload_evidence_file(mq_evidence_file, raw_file_name)

    assert data.shape == (58, 10), \
        "A number of extracted rows/columns is wrong."
    assert data['Raw file'].nunique() == 1 and \
        data['Raw file'].unique() == raw_file_name, \
        "Data not only for the specified raw file were extracted."
    assert sum(data['MS/MS scan number'].isna()) == 0, \
        "NA values in 'MS/MS scan number' column were not dropped."
