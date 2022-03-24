#!python
"""
This module provides functions to read MQ/DiaNN/AlphaPept output files and other IO supplementary functions.
"""

import logging
import os
import pandas as pd
import alphaviz.preprocessing


def read_file(
    filepath: str,
    column_names: list
) -> pd.DataFrame:
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
            if i == 0:  # for the first line to extract the index of the specified columns
                line = line.strip().split(sep)
                filename_col_index = [line.index(col) for col in column_names]
            else:  # use these indices for all other rows to extract the data
                line = line.split(sep)
                filename_data.append([line[ind] for ind in filename_col_index])
            i += 1

    data = pd.DataFrame(filename_data, columns=column_names)

    return data


def import_mq_evidence(
    filepath: str,
    experiment: str
) -> pd.DataFrame:
    """Read some columns from the output file evidence.txt of MaxQuant software.

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
            - 'Mass error [ppm]' ('float:.4d' type),
            - 'Modified sequence'.
        Renamed columns are marked as is the output data type of all columns. The rows of the data frame with missing 'MS/MS scan number' values are dropped.
    """
    chunk = pd.read_csv(filepath, chunksize=1000000, sep='\t', low_memory=False)
    data_raw_file = pd.concat(chunk)
    data_raw_file = data_raw_file[data_raw_file['Raw file'] == experiment]
    data_raw_file.rename(
        columns={
            'Score': 'Andromeda score',
            'K0': '1/K0',
        },
        inplace=True
    )
    data_raw_file.dropna(
        axis=0,
        subset=['MS/MS scan number', 'Proteins'],
        inplace=True
    )
    if 'Gene names' not in data_raw_file.columns:
        data_raw_file['Gene names'] = data_raw_file['Proteins'].apply(
            lambda x: ';'.join([entry.split('|')[-1].split('_')[0] for entry in x.split(';') if 'sp' in entry])
        )
    for col in ['Charge', 'MS/MS count', 'Gene names', 'Raw file']:
        data_raw_file[col] = data_raw_file[col].astype('category')
    for col in ['Retention time', 'Mass', 'm/z', '1/K0', 'Uncalibrated mass error [ppm]', 'Mass error [ppm]']:
        data_raw_file[col] = data_raw_file[col].astype(float).round(4)
    for col in ['MS/MS scan number', 'Length', 'Intensity', 'Andromeda score']:
        data_raw_file[col] = pd.to_numeric(
            data_raw_file[col],
            downcast='integer'
        )
    data_raw_file.dropna(
        axis=0,
        subset=['MS/MS scan number', 'Gene names'],
        inplace=True
    )
    first_column_names = ['Charge', 'm/z', 'Mass', '1/K0', 'Retention time']
    columns = list(data_raw_file.columns.drop(first_column_names))
    columns[1:1] = first_column_names
    data_raw_file = data_raw_file[columns]
    return data_raw_file


def import_mq_protein_groups(
    filepath: str,
    experiment: str
) -> pd.DataFrame:
    """Read the output file proteinGroups.txt of MaxQuant software.

    Parameters
    ----------
    filepath : str
        Full path to the proteinGroups.txt file.
    experiment : str
        The name of the experiment.

    Returns
    -------
    # pd.DataFrame
    #     The output data frame contains information about the following MQ columns:
    #         - 'Protein IDs',
    #         - 'Protein names',
    #         - 'Gene names',
    #         - 'Number of proteins' (renamed to '# proteins'),
    #         - 'Mol. weight [kDa]' (renamed to 'Mol weight, kDa'),
    #         - f'Peptides Exp_{experiment}' (renamed to '(EXP) # peptides'),
    #         - f'Unique peptides Exp_{experiment}' (renamed to '(EXP) # unique peptides'),
    #         - f'Sequence coverage Exp_{experiment} [%]' (renamed to '(EXP) Seq coverage, %'),
    #         - 'MS/MS count' (renamed to '# MS/MS'),
    #         - 'Sequence lengths',
    #     Renamed columns are marked. The rows of the data frame with missing 'Gene names' values are dropped.
    """
    data_common = pd.read_csv(filepath, sep='\t', low_memory=False)

    try:
        data_common.drop([col for col in data_common.columns if 'IDs' in col and col != 'Protein IDs'] + ['Best MS/MS', 'Peptide is razor'], inplace=True, axis=1)
    except:
        pass

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
        subset=['Fasta headers'],
        inplace=True
    )
    if '(EXP) # peptides' not in data_common.columns:
        data_common.rename(
            columns={
                'Peptides': '(EXP) # peptides',
            },
            inplace=True
        )
    try:
        data_common.dropna(
            axis=0,
            subset=['(EXP) # peptides'],
            inplace=True
        )
        data_common['(EXP) # peptides'] = data_common['(EXP) # peptides'].astype('int')
    except KeyError:
        pass

    try:
        data_common.Score = data_common.Score.astype(float)
    except:
        data_common.Score = data_common.Score.apply(lambda x: float(x) if x.replace('.', '', 1).isdigit() else None)

    if 'Gene names' not in data_common.columns:
        data_common[['Protein names', 'Protein IDs', 'Gene names']] = data_common.apply(lambda x: alphaviz.preprocessing.get_protein_info_from_fastaheader(x['Fasta headers']), axis=1, result_type='expand')
    data_common.dropna(
        axis=0,
        subset=['Gene names', 'Protein IDs', 'Score'],
        inplace=True
    )

    first_columns = [
        'Protein IDs',
        'Protein names',
        'Gene names',
        '# proteins',
        'Mol weight, kDa',
        '# MS/MS',
        'Sequence lengths',
    ]
    first_columns.extend([col for col in data_common.columns if '(EXP)' in col])

    data_common = data_common[first_columns + sorted(list(set(data_common.columns).difference(first_columns)))]

    return data_common


def import_mq_all_peptides(
    filepath: str
) -> pd.DataFrame:
    """Read some columns from the output file allPeptides.txt of MaxQuant software.

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


def import_mq_msms(
    filepath: str
) -> pd.DataFrame:
    """Read some columns from the output file msms.txt of MaxQuant software.

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
    try:
        maxquant_msms_columns = [
            'Scan number',
            'Matches',
            'Masses',
            'Mass deviations [Da]',
            'Mass deviations [ppm]'
        ]
        data_common = read_file(filepath, maxquant_msms_columns)
    except ValueError:
        maxquant_msms_columns = [
            'Scan number',
            'Matches',
            'Masses',
            'Mass Deviations [Da]',
            'Mass Deviations [ppm]'
        ]
        data_common = read_file(filepath, maxquant_msms_columns)
    data_common.columns = [col.strip().replace('Deviations', 'deviations') for col in data_common.columns]
    data_common['Scan number'] = data_common['Scan number'].astype('int')
    return data_common


def import_mq_summary(
    filepath: str
) -> pd.DataFrame:
    """Read the output file summary.txt of MaxQuant software.

    Parameters
    ----------
    filepath : str
        Full path to the msms.txt file.

    Returns
    -------
    pd.DataFrame
        The output data frame contains summary information of all the experiments.
    """
    data_common = pd.read_csv(filepath, sep='\t', low_memory=False)
    data_common.dropna(subset=['MS'], axis=0, inplace=True)
    return data_common


def import_mq_output(
    necessary_files: list,
    path_mq_output_folder: str,
    experiment: str
):
    """Read all specified files from the MQ output folder and returns the data frames for each of the files.

    Parameters
    ----------
    necessary_files : list
        A list of strings containing the names of the MQ output files with extensions, e.g. ['allPeptides.txt', 'msms.txt'].
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
        'allPeptides.txt': import_mq_all_peptides,
        'msms.txt': import_mq_msms,
        'evidence.txt': import_mq_evidence,
        'proteinGroups.txt': import_mq_protein_groups,
        'summary.txt': import_mq_summary,
    }
    for file in necessary_files:
        file_path = os.path.join(
            path_mq_output_folder,
            file
        )
        if file in ['allPeptides.txt', 'msms.txt', 'summary.txt']:
            df = file_func_dict[file](
                file_path
            )
        else:
            df = file_func_dict[file](
                file_path,
                experiment
            )
        logging.info(f"MaxQuant output {file} file is uploaded.")
        yield df


def get_filenames_from_directory(
    directory: str,
    extensions_list: list
) -> list:
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


def read_fasta(
    filepath: str
) -> dict:
    """Read the fasta file using the pyteomics package.

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


def import_diann_stats(
    filepath: str,
    experiment: str
):
    """Read the DIANN output .stats.tsv file.

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
    diann_overview = pd.read_csv(filepath, sep='\t', low_memory=False)
    return diann_overview


def create_diann_proteins_table(
    diann_df: pd.DataFrame,
    fasta: object
):
    """Extract information about genes, proteins and protein groups from the loaded main DIANN output .tsv file.

    Parameters
    ----------
    diann_df : pd.DataFrame
        The original data frame after loading the main .tsv DIANN output file and filter by the experiment name.
    fasta : pyteomics.fasta.IndexedUniProt
        The object containing information about all proteins from the fasta file.

    Returns
    -------
    pd.DataFrame
        The output data frame contains information about genes, proteins and proteins groups.
    """
    columns = [col for col in diann_df.columns if 'PG' in col or 'Protein' in col or 'Genes' in col]
    cols_to_remove = ['Protein.Group', 'Protein.Ids', 'Protein.Names']
    for col in cols_to_remove:
        columns.remove(col)
    proteins = diann_df.groupby(columns).agg({
        'Protein.Ids': lambda x: ','.join(set(x)),
        'MS2.Scan': lambda x: len(set(x)),
        'Stripped.Sequence': lambda x: len(set(x))
    }).reset_index()
    proteins.rename(columns={
        'MS2.Scan': '# MS/MS',
        'Stripped.Sequence': '(EXP) # peptides',
        'Genes': 'Gene names',
        'Protein.Ids': 'Protein IDs'
    }, inplace=True)
    proteins['# proteins'] = proteins['Protein IDs'].apply(lambda x: len(x.split(',')))
    proteins['Protein names'], proteins['Sequence lengths'] = zip(
        *proteins['Protein IDs'].apply(lambda x: alphaviz.preprocessing.get_protein_info(fasta, x)))
    first_columns = ['Protein IDs', 'Protein names', 'Gene names', '# proteins', '(EXP) # peptides', '# MS/MS', 'Sequence lengths']
    proteins = proteins[first_columns + sorted(list(set(proteins.columns).difference(first_columns)))]
    return proteins


def create_diann_peptides_table(
    diann_df: pd.DataFrame
):
    """Extract information about peptides from the loaded main DIANN output .tsv file.

    Parameters
    ----------
    diann_df : pd.DataFrame
        The original data frame after loading the main .tsv DIANN output file and filter by the experiment name.

    Returns
    -------
    pd.DataFrame
        The output data frame contains information about peptides.
    """
    peptides = diann_df.copy()
    columns = [col for col in peptides.columns if 'PG' not in col and 'Protein' not in col and 'Genes' not in col and 'GG' not in col]
    columns.extend(['Genes'])

    peptides = diann_df[columns[2:]].copy()
    peptides['Length'] = peptides['Stripped.Sequence'].str.len()

    peptides.rename(columns={
        'MS2.Scan': 'MS/MS scan number',
        'Genes': 'Gene names',
        'Precursor.Charge': 'Charge',
        'Stripped.Sequence': 'Sequence'
    }, inplace=True)

    peptides['Sequence_AP_mod'] = peptides['Modified.Sequence'].apply(alphaviz.preprocessing.convert_diann_ap_mod)
    peptides['Modified.Sequence'] = peptides['Modified.Sequence'].apply(alphaviz.preprocessing.convert_diann_mq_mod)
    peptides['m/z'] = 0.0
    first_columns = ['Modified.Sequence', 'Length', 'm/z', 'RT', 'Predicted.RT', 'Charge', 'IM', 'Predicted.IM']
    peptides = peptides[first_columns + sorted(list(set(peptides.columns).difference(first_columns)))]
    return peptides


def import_diann_output(
    path_diann_output_folder: str,
    experiment: str,
    fasta: object
):
    """Load two files from the DiaNN output folder and returns the data frames containing information about proteins, peptides, and summary information about the whole experiment.

    Parameters
    ----------
    path_diann_output_folder : str
        Path to the DIANN output folder with all output files needed.
    experiment : str
        The name of the experiment.
    fasta : pyteomics.fasta.IndexedUniProt
        The object containing information about all proteins from the fasta file.

    Returns
    -------
    list of pd.DataFrames
        The function returns three pandas data frame with the extracted information about proteins, peptides, and summary information about the whole experiment.
    """
    diann_output_file, diann_stats_file = sorted(get_filenames_from_directory(
        path_diann_output_folder, 'tsv'), key=len)[:2]

    diann_df = pd.read_csv(os.path.join(path_diann_output_folder, diann_output_file), sep='\t', low_memory=False)
    diann_df = diann_df[diann_df.Run == experiment]

    diann_proteins = create_diann_proteins_table(diann_df, fasta)
    diann_peptides = create_diann_peptides_table(diann_df)

    diann_overview = import_diann_stats(os.path.join(path_diann_output_folder, diann_stats_file), experiment)

    return diann_proteins, diann_peptides, diann_overview, diann_output_file
