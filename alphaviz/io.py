#!python
"""
This module provides functions to read MQ/DiaNN/AlphaPept output files and other IO supplementary functions.
"""

import logging
import os
import pandas as pd


def read_file(
    filepath: str,
    column_names: list
)-> pd.DataFrame:
    """Enable reading the file and retrieving the values from the
    specified columns. Compared to function pd.read_csv() it gains significant time if the file is huge and is only a few ms slower for small files.

    Parameters
    ----------
    filepath : str
        Full path to the file.
    column_names : list
        A list of column names to be read.

    Returns
    -------
    pd.DataFrame
        This data frame contains data from all columns of the specified file.
    """
    file_ext = os.path.splitext(filepath)[-1]
    if file_ext == '.csv':
        sep = ','
    elif file_ext in ['.tsv', '.txt']:
        sep = '\t'
    with open(filepath) as filelines:
        i = 0
        filename_col_index = []
        filename_data = []

        for line in filelines:
            if i == 0: # for the first line to extract the index of the specified columns
                line = line.strip().split(sep)
                filename_col_index = [line.index(col) for col in column_names]
            else: # use these indices for all other rows to extract the data
                line = line.split(sep)
                filename_data.append([line[ind] for ind in filename_col_index])
            i += 1

    data = pd.DataFrame(filename_data, columns=column_names)

    return data


def upload_evidence_file(
    filepath: str,
    experiment: str
)-> pd.DataFrame:
    """Load some columns from the output file evidence.txt of MaxQuant software.

    Parameters
    ----------
    filepath : str
        Full path to the evidence.txt file.
    experiment : str
        The name of the experiment.

    Returns
    -------
    pd.DataFrame
        The output data frame contains information about the following MQ columns:
            - 'Sequence',
            - 'Length' ('int' type),
            - 'Acetyl (Protein N-term)' (renamed to 'Acetylation (N-term)') ('int' type),
            - 'Oxidation (M)' ('int' type),
            - 'Proteins',
            - 'Retention time' ('float:.4d' type),
            - 'Mass' ('float:.4d' type),
            - 'm/z' ('float:.4d' type),
            - 'Charge' ('category' type),
            - 'Intensity' ('int' type),
            - '1/K0' ('float:.4d' type),
            - 'MS/MS count' ('category' type),
            - 'MS/MS scan number' ('int' type),
            - 'Gene names' ('category' type),
            - 'Score' (renamed to 'Andromeda score') ('int' type),
            - 'Raw file' ('category' type),
            - 'Uncalibrated mass error [ppm]' ('float:.4d' type),
            - 'Mass error [ppm]' ('float:.4d' type).
        Renamed columns are marked as is the output data type of all columns. The rows of the data frame with missing 'MS/MS scan number' values are dropped.
    """
    maxquant_evidence_columns = [
        'Sequence',
        'Length',
        'Acetyl (Protein N-term)',
        'Oxidation (M)',
        'Proteins',
        'Retention time',
        'Mass',
        'm/z',
        'Charge',
        'Intensity',
        '1/K0',
        'MS/MS count',
        'MS/MS scan number',
        'Gene names',
        'Score',
        'Raw file',
        'Uncalibrated mass error [ppm]',
        'Mass error [ppm]'
    ]
    data_common = read_file(filepath, maxquant_evidence_columns)
    data_common.rename(
        columns={
            'Acetyl (Protein N-term)': 'Acetylation (N-term)',
            'Score': 'Andromeda score'
        },
        inplace=True
    )
    data_raw_file = data_common[data_common['Raw file'] == experiment]

    for col in ['Charge', 'MS/MS count', 'Gene names', 'Raw file']:
        data_raw_file[col] = data_raw_file[col].astype('category')
    for col in ['Retention time', 'Mass', 'm/z', '1/K0', 'Uncalibrated mass error [ppm]', 'Mass error [ppm]']:
        data_raw_file[col] = data_raw_file[col].astype(float).round(4)
    for col in ['MS/MS scan number', 'Acetylation (N-term)', 'Oxidation (M)', 'Length', 'Intensity', 'Andromeda score']:
        data_raw_file[col] = pd.to_numeric(
            data_raw_file[col],
            downcast='integer'
        )
    data_raw_file.dropna(
        axis=0,
        subset=['MS/MS scan number']
    )
    return data_raw_file


def upload_protein_groups_file(
    filepath: str,
    experiment: str
)-> pd.DataFrame:
    """Load some columns from the output file proteinGroups.txt of MaxQuant software.

    Parameters
    ----------
    filepath : str
        Full path to the proteinGroups.txt file.
    experiment : str
        The name of the experiment.

    Returns
    -------
    pd.DataFrame
        The output data frame contains information about the following MQ columns:
            - 'Protein IDs',
            - 'Protein names',
            - 'Gene names',
            - 'Number of proteins' (renamed to '# proteins'),
            - 'Mol. weight [kDa]' (renamed to 'Mol weight, kDa'),
            - f'Peptides Exp_{experiment}' (renamed to '(EXP) # peptides'),
            - f'Unique peptides Exp_{experiment}' (renamed to '(EXP) # unique peptides'),
            - f'Sequence coverage Exp_{experiment} [%]' (renamed to '(EXP) Seq coverage, %'),
            - 'MS/MS count' (renamed to '# MS/MS'),
            - 'Sequence lengths',
        Renamed columns are marked. The rows of the data frame with missing 'Gene names' values are dropped.
    """
    maxquant_protein_groups_columns = [
        'Protein IDs',
        'Protein names',
        'Gene names',
        'Number of proteins',
        'Mol. weight [kDa]',
        f'Peptides Exp_{experiment}',
        f'Unique peptides Exp_{experiment}',
        f'Sequence coverage Exp_{experiment} [%]',
        'MS/MS count',
        'Sequence lengths',
    ]
    data_common = read_file(filepath, maxquant_protein_groups_columns)
    data_common.rename(
        columns={
            'Number of proteins': '# proteins',
            'Mol. weight [kDa]': 'Mol weight, kDa',
            f'Peptides Exp_{experiment}': '(EXP) # peptides',
            f'Unique peptides Exp_{experiment}': '(EXP) # unique peptides',
            f'Sequence coverage Exp_{experiment} [%]': '(EXP) Seq coverage, %',
            'MS/MS count': '# MS/MS'
        },
        inplace=True
    )
    data_common.dropna(
        axis=0,
        subset=['Gene names'],
        inplace=True
    )
    return data_common


def upload_all_peptides_file(
    filepath:str
)-> pd.DataFrame:
    """Load some columns from the output file allPeptides.txt of MaxQuant software.

    Parameters
    ----------
    filepath : str
        Full path to the allPeptides.txt file.

    Returns
    -------
    pd.DataFrame
        The output data frame contains information about the following MQ columns:
            - 'Pasef MS/MS IDs' ('list' type),
            - 'MS/MS scan number' ('int' type).
        The rows of the data frame with missing 'MS/MS scan number' values are dropped.
    """
    maxquant_all_peptides_columns = [
        'Pasef MS/MS IDs',
        'MS/MS scan number'
    ]
    data_common = read_file(filepath, maxquant_all_peptides_columns)
    data_common.columns = [col.strip() for col in data_common.columns]
    data_common['MS/MS scan number'] = data_common['MS/MS scan number'].str.strip()
    data_common = data_common[data_common['MS/MS scan number'] != '']
    data_common['MS/MS scan number'] = data_common['MS/MS scan number'].astype(int)
    data_common['Pasef MS/MS IDs'] = data_common['Pasef MS/MS IDs'].str.split(';')
    return data_common


def upload_msms_file(
    filepath: str
)-> pd.DataFrame:
    """Load some columns from the output file msms.txt of MaxQuant software.

    Parameters
    ----------
    filepath : str
        Full path to the msms.txt file.

    Returns
    -------
    pd.DataFrame
        The output data frame contains information about the following MQ columns:
            - 'Scan number' ('int' type),
            - 'Matches',
            - 'Masses',
            - 'Mass deviations [Da]',
            - 'Mass deviations [ppm]'.
    """
    maxquant_msms_columns = [
        'Scan number',
        'Matches',
        'Masses',
        'Mass deviations [Da]',
        'Mass deviations [ppm]'
    ]
    data_common = read_file(filepath, maxquant_msms_columns)
    data_common.columns = [col.strip() for col in data_common.columns]
    data_common['Scan number'] = data_common['Scan number'].astype('int')
    return data_common


def upload_mq_files(
    necessary_files: list,
    path_mq_output_folder: str,
    experiment: str
):
    """Load all specified files from the MQ output folder and returns the data frames for each of the files.

    Parameters
    ----------
    necessary_files : list
        A list of strings containing the names of the MQ output files without extensions, e.g. ['allPeptides', 'msms'].
    path_mq_output_folder : str
        Path to the MaxQuant output folder with all output files needed.
    experiment : str
        The name of the experiment.

    Returns
    -------
    generator
        For each of the specified MQ output files, the function returns a pandas data frame with the extracted information.
    """
    file_func_dict = {
        'allPeptides': upload_all_peptides_file,
        'msms': upload_msms_file,
        'evidence': upload_evidence_file,
        'proteinGroups': upload_protein_groups_file
    }
    for file in necessary_files:
        file_path = os.path.join(
            path_mq_output_folder,
            f'{file}.txt'
        )
        if file in ['allPeptides', 'msms']:
            df = file_func_dict[file](
                file_path
            )
        else:
            df = file_func_dict[file](
                file_path,
                experiment
            )
        logging.info(f"MaxQuant output {file}.txt file is uploaded.")
        yield df


def get_file_names_from_directory(
    directory: str,
    extensions_list: list
)-> list:
    """Search for files with the specified extension in the repository and return a list of all file names with that extention.

    Parameters
    ----------
    directory : str
        Path to the repository to search in.
    extensions_list : list
        A list of extensions, e.g. ['d', 'hdf'].

    Returns
    -------
    list
        The list of filtered file names based on their extensions.
    """
    file_names = [file for file in os.listdir(directory) if file.split('.')[-1] in extensions_list]
    return file_names


def upload_fasta_file(
    filepath: str
)-> dict:
    """Load the fasta file using the pyteomics package.

    Parameters
    ----------
    filepath : str
        Full path to the .fasta file.

    Returns
    -------
    pyteomics.fasta.IndexedUniProt object
        The output object allows access to all available information in the fasta file using the protein ID.
    """
    import pyteomics.fasta
    fasta = pyteomics.fasta.IndexedUniProt(filepath)
    return fasta


def load_diann_stats_file(
    filepath: str,
    experiment: str
):
    """Load the DIANN output .stats.tsv file.

    Parameters
    ----------
    filepath : str
        Full path to the .stats.tsv file.
    experiment : str
        The name of the experiment.

    Returns
    -------
    pd.DataFrame
        The output data frame contains summary information about the whole experiment.
    """
    diann_overview = pd.read_csv(filepath, sep='\t')
    diann_overview = diann_overview[diann_overview['File.Name'].str.contains(experiment)]
    diann_overview = diann_overview[diann_overview.columns[1:]].T
    diann_overview.reset_index(inplace=True)
    diann_overview.rename(columns={'index': 'parameters', 0: 'values'}, inplace=True)
    diann_overview['values'] = diann_overview['values'].apply(lambda x: '%.2E' % x if x>100000 else '%.2f' % x)
    return diann_overview
