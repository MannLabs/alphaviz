#!python
"""
This module provides functions that are helping to preprocess the data.
"""

import re
import logging
import pandas as pd


def get_mq_unique_proteins(
    filepath: str
) -> list:
    """Extract unique "Protein names" from the specified MaxQuant output file.

    Parameters
    ----------
    filepath : str
        Full path to the file.

    Returns
    -------
    list
        A list of unique protein names from the specified file.

    """
    if filepath.endswith('.txt'):
        sep = '\t'
    else:
        raise ValueError('The given file does not have an original .txt extension.')

    with open(filepath) as filelines:
        i = 0
        filename_col_index = int()
        filename_data = []

        for line in filelines:
            line = line.split(sep)
            if i == 0:
                try:
                    filename_col_index = line.index('Protein names')
                except ValueError:
                    raise ValueError('The given file does not contain the "Protein names" column.')
            else:
                filename_data.append(line[filename_col_index])
            i += 1

        unique_proteins = set(filename_data)

    sorted_unique_proteins = sorted(list(unique_proteins))
    return sorted_unique_proteins


def filter_df(
    df: pd.DataFrame,
    pattern: str,
    column: str,
    software: str
) -> pd.DataFrame:
    """Filter the data frame based on the pattern (any value) in the specified column.

    Parameters
    ----------
    df : pd.DataFrame
        The original data frame.
    pattern : str
        The string to be used to filter of the data frame column.
    column : str
        The column to be used to filter.
    software: str
        The name of the software tool where the filtering is used.

    Returns
    -------
    pd.DataFrame
        The filtered data frame.

    """
    if not pattern:
        return df
    if software == 'maxquant':
        output = df[df[column].str.contains(pattern, na=False)]
    else:
        output = df[df[column] == pattern]
    return output


def sort_naturally(
    line: str,
    reverse: bool = False
) -> str:
    """Sort the string natural to humans, e.g. 4,1,6,11 will be sorted as 1,4,6,11 and not like 1,11,4,6.

    Parameters
    ----------
    line : str
        The string to be sorted.
    reverse : bool
        Whether to apply the reverse option or not. Defaults: False.

    Returns
    -------
    str
        A naturally sorted string.

    """
    def _convert(value):
        if value.isdigit():
            return int(value)
        return value.lower()
    alphanum_key = lambda key: [_convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(line, key=alphanum_key, reverse=reverse)


def get_aa_seq(
    protein_id: str,
    fasta,  # pyteomics.fasta.IndexedUniProt object
) -> str:
    """Extract the leading razor protein sequence for the list of
    proteinIDs of the protein group from the pyteomics.fasta.IndexedUniProt object.

    Parameters
    ----------
    protein_ids : str
        String containing the proteinIDs of all protein isoforms.
    fasta : pyteomics.fasta.IndexedUniProt object
        The pyteomics.fasta.IndexedUniProt object.

    Returns
    -------
    str
        Protein sequence for the leading razor protein, e.g. from the list of proteinIDs 'Q15149;Q15149-7;Q15149-9;Q15149-5' the AA sequence for protein Q15149 will be returned.

    """
    # for id in sorted(protein_ids.split(';'), reverse=True):
    try:
        protein_seq = fasta.get_by_id(protein_id).sequence
        return protein_seq
    except KeyError:
        logging.info(f"The provided protein ID {id} is missing in the fasta file.")


def get_mq_ms2_scan_data(
    msms: pd.DataFrame,
    selected_msms_scan: int,
    raw_data,  # AlphaTims TimsTOF object,
    precursor_id: int
) -> pd.DataFrame:
    """Extract MS2 data as a data frame for the specified MSMS scan number and precursor ID from the 'msms.txt' MQ output file and raw file.

    Parameters
    ----------
    msms : pd.DataFrame
        Pre-loaded 'msms.txt' MQ output file.
    selected_msms_scan : int
        MSMS scan number.
    raw_data : AlphaTims TimsTOF object
        AlphaTims TimsTOF object.
    precursor_id : int
        The identifier of the precursor.

    Returns
    -------
    pd.DataFrame
        For the specified MSMS scan and precursor ID, the extracted data frame contains the following columns:
            - 'mz_values'
            - 'intensity_values'
            - 'ions'
            - 'wrong_dev_value': whether the mass_deviation specified in the MQ table was incorrect.

    """
    msms_filtered = msms[msms['Scan number'] == selected_msms_scan]
    msms_filtered_df = pd.DataFrame(
        {
            'ions': msms_filtered.Matches.values[0].split(';'),
            'mz': msms_filtered.Masses.values[0].split(';'),
            'mass_dev_Da': msms_filtered['Mass deviations [Da]'].values[0].split(';'),
            'mass_dev_ppm': msms_filtered['Mass deviations [ppm]'].values[0].split(';')
        }
    )
    for col in ['mz', 'mass_dev_Da', 'mass_dev_ppm']:
        msms_filtered_df[col] = msms_filtered_df[col].astype(float)

    data = raw_data[:, :, precursor_id].loc[:, ['mz_values', 'intensity_values']]  # can be slightly faster by only retrieving the indices and converting directly to mz values and intensities
    data['ions'] = '-'
    data['wrong_dev_value'] = False

    for row in msms_filtered_df.itertuples():  # inefficient implementation. Sorting both arrays by mz should allow you to do it much faster.
        ion_index = data.mz_values.sub(row.mz).abs().idxmin()
        mass_dev_ppm_calc = ((row.mz + row.mass_dev_Da - data.loc[ion_index, 'mz_values']) * 10**6) / data.loc[ion_index, 'mz_values']
        # to think how not to set a fixed threshold?
        ppm_threshold = 100
        if abs(mass_dev_ppm_calc) < ppm_threshold:
            data.loc[ion_index, 'ions'] = row.ions
            msms_filtered_df.loc[msms_filtered_df.ions == row.ions, 'mass_dev_ppm'] = mass_dev_ppm_calc

    data.drop_duplicates('mz_values', inplace=True)
    data.sort_values(['ions', 'intensity_values'], ascending=True, inplace=True)
    data_merged = pd.merge(data, msms_filtered_df, on='ions', how='left')

    return data_merged.drop('mz', axis=1)


def get_identified_ions(
    values: list,
    sequence: str,
    ion_type: str
) -> list:
    """For the specified peptide sequence extract all identified in the experiment ions and based on the specified ion_type return a list of booleans containing information for the b-ions whether the peptide is breaking after aligned amino acid or for the y-ion whether is breaking before aligned amino acid.

    E.g. for the peptide 'NTINHN' having the unique ion values ['b2-H2O', 'b2', 'b3', 'b5-NH3'] it will return the following list of booleans: [False,True,True,False,True,False].

    Parameters
    ----------
    values : list
        The list of all unique identified ions for the peptide in the experiment.
    sequence : str
        Peptide sequence.
    ion : str
        Ion type, e.g. 'b' or 'y'. Other ion types are not implemented.

    Returns
    -------
    list
        List of peptide length of booleans with True for the presenting ion and False for a missing one.

    """
    ions = [False] * (len(sequence))
    all_ions = list(set([ion.split('-')[0] for ion in set(values) if ion_type in ion]))
    for each in all_ions:
        if ion_type == 'b':
            ions[int(each.replace(ion_type, '')) - 1] = True
        elif ion_type == 'y':
            ions[-int(each.replace(ion_type, ''))] = True
        else:
            raise NotImplementedError(f"The specified ion type {ion_type} is not implemented in the current version.")
    return ions


def convert_diann_mq_mod(
    sequence: str
) -> str:
    # this function is taken from the AlphaMap package and modified
    """Convert DIA-NN style modifications to MaxQuant style modifications.

    Args:
        sequence (str): A peptide sequence with a DIA-NN style modification.

    Returns:
        str: A peptide sequence with MaxQuant style modification.
    """

    modif_convers_dict = {
        '(UniMod:1)': '[Acetyl ({})]',
        '(UniMod:2)': '[Amidated ({})]',
        '(UniMod:4)': '[Carbamidomethyl ({})]',
        '(UniMod:5)': '[Carbamyl ({})]',
        '(UniMod:7)': '[Deamidation ({})]',
        '(UniMod:21)': '[Phospho ({})]',
        '(UniMod:23)': '[Dehydrated ({})]',
        '(UniMod:26)': '[Pyro-carbamidomethyl ({})]',
        '(UniMod:27)': '[Glu->pyro-Glu]',
        '(UniMod:28)': '[Gln->pyro-Glu]',
        '(UniMod:30)': '[Cation:Na ({})]',
        '(UniMod:34)': '[Methyl ({})]',
        '(UniMod:35)': '[Oxidation ({})]',
        '(UniMod:36)': '[Dimethyl ({})]',
        '(UniMod:37)': '[Trimethyl ({})]',
        '(UniMod:40)': '[Sulfo ({})]',
        '(UniMod:55)': '[Cys-Cys]',
        '(UniMod:121)': '[GlyGly ({})]',
        '(UniMod:254)': '[Delta:H(2)C(2) ({})]',
        '(UniMod:312)': '[Cysteinyl]',
        '(UniMod:345)': '[Trioxidation ({})]',
        '(UniMod:408)': '[Hydroxyproline]',
        '(UniMod:425)': '[Dioxidation ({})]',
        '(UniMod:526)': '[Dethiomethyl ({})]',
        '(UniMod:877)': '[QQTGG ({})]',
    }
    mods = re.findall(r'\(UniMod:\d+\)', sequence)
    if mods:
        for mod in mods:
            posit = re.search(r'\(UniMod:\d+\)', sequence)
            i = posit.start()

            if i == 0:
                add_aa = 'N-term'
            elif posit.end() == len(sequence):
                add_aa = 'C-term'
            else:
                add_aa = sequence[i-1]

            if mod == '(UniMod:7)':
                if add_aa in 'NQ':
                    add_aa = 'NQ'
            elif mod == '(UniMod:21)':
                if add_aa in 'STY':
                    add_aa = 'STY'
            elif mod == '(UniMod:23)':
                if add_aa in 'ST':
                    add_aa = 'ST'
            elif mod == '(UniMod:30)':
                if add_aa in 'DE':
                    add_aa = 'DE'
            elif mod == '(UniMod:34)':
                if add_aa in 'KR':
                    add_aa = 'KR'
            elif mod == '(UniMod:36)':
                if add_aa in 'KR':
                    add_aa = 'KR'
            elif mod == '(UniMod:40)':
                if add_aa in 'STY':
                    add_aa = 'STY'
            elif mod == '(UniMod:425)':
                if add_aa in 'MW':
                    add_aa = 'MW'

            if mod in modif_convers_dict.keys():
                sequence = sequence.replace(mod, modif_convers_dict.get(mod).format(add_aa), 1)
            else:
                logging.info(f"This modification {mod} can't be converted.")

    return sequence


def convert_diann_ap_mod(
    sequence: str
) -> str:
    """Convert DIA-NN style modifications to AlphaPept style modifications.

    Args:
        sequence (str): A peptide sequence with a DIA-NN style modification.

    Returns:
        str: A peptide sequence with AlphaPept style modification.
    """
    modif_convers_dict = {
        '(UniMod:1)': 'a',  # '[Acetyl ({})]'
        '(UniMod:2)': 'am',  # '[Amidated ({})]'
        '(UniMod:4)': 'c',  # '[Carbamidomethyl ({})]'
        '(UniMod:7)': 'deam',  # '[Deamidation ({})]'
        '(UniMod:21)': 'p',  # '[Phospho ({})]'
        '(UniMod:26)': 'cm',  # '[Pyro-carbamidomethyl ({})]',
        '(UniMod:27)': 'pg',  # '[Glu->pyro-Glu]'
        '(UniMod:28)': 'pg',  # '[Gln->pyro-Glu]'
        '(UniMod:35)': 'ox',  # '[Oxidation ({})]'
    }
    mods = re.findall(r'\(UniMod:\d+\)', sequence)

    if mods:
        for mod in mods:
            posit = re.search(r'\(UniMod:\d+\)', sequence)
            i = posit.start()
            if i != 0:
                i -= 1
            if mod in modif_convers_dict.keys():
                sequence = sequence.replace(mod, '', 1)
                sequence = sequence[:i] + modif_convers_dict[mod] + sequence[i:]
            else:
                logging.info(f"This modification {mod} can't be converted.")

    return sequence


def get_protein_info(
    fasta: dict,
    protein_ids: str
):
    """Get the name and the length of the protein(s) from the fasta file specifying the protein id(s).

    Parameters
    ----------
    fasta : pyteomics.fasta.IndexedUniProt object
        The Pyteomics object contains information about all proteins from the .fasta file.
    protein_ids : str
        The list of the protein IDs separated by comma.

    Returns
    -------
    tuple of strings
        The name and the length of the specified protein(s).

    """
    protein_names = []
    protein_seq_lens = []
    for protein_id in protein_ids.split():
        try:
            protein_names.append(fasta.get_by_id(protein_id).description['name'])
        except KeyError:
            logging.info(f"The protein id {protein_id} is not found in the fasta file.")
        try:
            protein_seq_lens.append(str(len(fasta.get_by_id(protein_id).sequence)))
        except KeyError:
            logging.info(f"The sequence length for the protein {protein_id} is not found in the fasta file.")
    return ','.join(protein_names), ','.join(protein_seq_lens)


def get_protein_info_from_fastaheader(
    string: str,
    **kwargs
):
    """Extract information about protein IDs, protein names and gene names from the "Fasta headers" column of the MQ output tables.

    Parameters
    ----------
    string : str
        A 'Fasta header' string from the MQ output table for one protein group (e.g. from the proteinGroups.txt file).

        E.g. a complex one: 'sp|Q3SY84|K2C71_HUMAN Keratin, type II cytoskeletal 71 OS=Homo sapiens OX=9606 GN=KRT71 PE=1 SV=3;;sp|Q14CN4|K2C72_HUMAN Keratin, type II cytoskeletal 72 OS=Homo sapiens OX=9606 GN=KRT72 PE=1 SV=2;;;sp|Q7RTS7|K2C74_HUMAN Keratin, type II cytoskeletal 74 OS'

    Returns
    -------
    a tuple of strings
        The function returns a tuple of three strings containing information about the protein names, protein IDs and gene names.

    """
    protein_names = []
    protein_ids = []
    genes = []
    for protein in string.split(';'):
        if protein:
            try:
                protein_names.append(re.findall(r'\s(.+)OS', protein)[0].strip())
            except:
                protein_names.append("")
            try:
                protein_ids.append(protein.split('|')[1])
            except:
                protein_ids.append("")
            try:
                genes.append(protein.split()[0].split('|')[-1].split('_')[0])
            except:
                genes.append("")
    protein_names = ";".join(protein_names) if protein_names else None
    protein_ids = ";".join(protein_ids) if protein_ids else None
    genes = ";".join(genes) if genes else None
    return protein_names, protein_ids, genes
