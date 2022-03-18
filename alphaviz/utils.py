import os
import numba
from numba import types
from numba.typed import Dict, List
from numba import njit
import numpy as np

# paths
BASE_PATH = os.path.dirname(__file__)
IMG_PATH = os.path.join(BASE_PATH, "img")
STYLE_PATH = os.path.join(BASE_PATH, "style")
DOCS_PATH = os.path.join(BASE_PATH, "docs")
DATA_PATH = os.path.join(BASE_PATH, "data")
MODELS_PATH = os.path.join(BASE_PATH, "models")
LATEST_GITHUB_INIT_FILE = "https://github.com/MannLabs/alphaviz/blob/main/alphaviz/__init__.py"


def check_analysis_file(file):
    # TODO: write the checks for the preloaded file and return the exception when the file can't be uploaded
    pass


# this code was taken from the AlphaTims Python package (https://github.com/MannLabs/alphatims/blob/master/alphatims/utils.py) and modified
def check_github_version(silent=False) -> str:
    """Checks and logs the current version of AlphaViz.
    Check if the local version equals the AlphaViz GitHub master branch.
    This is only possible with an active internet connection and
    if no credentials are required for GitHub.
    Parameters
    ----------
    silent : str
        Use the logger to display the obtained conclusion.
        Default is False.
    Returns
    -------
    : str
        The version on the AlphaViz GitHub master branch.
        "" if no version can be found on GitHub
    """
    import requests
    from bs4 import BeautifulSoup
    import alphaviz

    try:
        main_response = requests.get(LATEST_GITHUB_INIT_FILE)
        main_soap = BeautifulSoup(main_response.content.decode('utf-8'), 'html.parser')
        for line in main_soap.find_all('td', class_='blob-code blob-code-inner js-file-line'):
            if line.text.startswith('__version__'):
                github_version = line.text.split()[-1].strip()[1:-1]
                if not silent:
                    if github_version != alphaviz.__version__:
                        print(f"You are currently using AlphaViz version {alphaviz.__version__}. However, the latest version of AlphaViz on GitHub is {github_version}. Checkout https://github.com/MannLabs/alphaviz.git for instructions on how to update AlphaViz...")
                    else:
                        print("Current AlphaViz version is up-to-date with GitHub.")
                return github_version
    except:
        print("Could not check GitHub for the latest AlphaViz release.")
        return ""


# This code was taken from the AlphaPept Python package (https://github.com/MannLabs/alphapept/blob/master/nbs/03_fasta.ipynb)
def get_mass_dict(
    modfile: str = "data/modifications.tsv",
    aasfile: str = "data/amino_acids.tsv",
    verbose: bool = True
):
    """
    Function to create a mass dict based on tsv files.
    This is used to create the hardcoded dict in the constants notebook.
    The dict needs to be hardcoded because of importing restrictions when using numba.
    More specifically, a global needs to be typed at runtime.
    Args:
        modfile (str): Filename of modifications file.
        aasfile (str): Filename of AAs file.
        verbose (bool, optional): Flag to print dict.
    Returns:
        Returns a numba compatible dictionary with masses.
    Raises:
        FileNotFoundError: If files are not found.
    """
    import pandas as pd

    mods = pd.read_csv(modfile, delimiter="\t")
    aas = pd.read_csv(aasfile, delimiter="\t")

    mass_dict = Dict.empty(key_type=types.unicode_type, value_type=types.float64)

    for identifier, mass in aas[["Identifier", "Monoisotopic Mass (Da)"]].values:
        mass_dict[identifier] = float(mass)

    for identifier, aar, mass in mods[
        ["Identifier", "Amino Acid Residue", "Monoisotopic Mass Shift (Da)"]
    ].values:

        if ("<" in identifier) or (">" in identifier):
            for aa_identifier, aa_mass in aas[["Identifier", "Monoisotopic Mass (Da)"]].values:
                if '^' in identifier:
                    new_identifier = identifier[:-2] + aa_identifier
                    mass_dict[new_identifier] = float(mass) + mass_dict[aa_identifier]
                elif aar == aa_identifier:
                    new_identifier = identifier[:-2] + aa_identifier
                    mass_dict[new_identifier] = float(mass) + mass_dict[aa_identifier]
                else:
                    pass
        else:
            mass_dict[identifier] = float(mass) + mass_dict[aar]

    # Manually add other masses
    mass_dict["Electron"] = (0.000548579909070)  # electron mass, half a millimass error if not taken into account
    mass_dict["Proton"] = 1.00727646687  # proton mass
    mass_dict["Hydrogen"] = 1.00782503223  # hydrogen mass
    mass_dict["C13"] = 13.003354835  # C13 mass
    mass_dict["Oxygen"] = 15.994914619  # oxygen mass
    mass_dict["OH"] = mass_dict["Oxygen"] + mass_dict["Hydrogen"]  # OH mass
    mass_dict["H2O"] = mass_dict["Oxygen"] + 2 * mass_dict["Hydrogen"]  # H2O mass

    mass_dict["NH3"] = 17.03052
    mass_dict["delta_M"] = 1.00286864
    mass_dict["delta_S"] = 0.0109135
    mass_dict['a'] = 42.01056469

    if verbose:

        for element in mass_dict:
            print('mass_dict["{}"] = {}'.format(element, mass_dict[element]))

    return mass_dict


@njit
def parse(
    peptide: str
) -> List:
    """
    Parser to parse peptide strings
    Args:
        peptide (str): modified peptide sequence.
    Return:
        List (numba.typed.List): a list of animo acids and modified amono acids
    """
    if "_" in peptide:
        peptide = peptide.split("_")[0]
    parsed = List()
    string = ""

    for ind, i in enumerate(peptide):
        string += i
        if ind == 0 and i == 'a' and peptide[1].islower():  # protein N-term modification
            parsed.append(string)
            string = ""
        if i.isupper():
            parsed.append(string)
            string = ""

    return parsed


@njit
def get_precmass(
    parsed_pep: list,
    mass_dict: numba.typed.Dict
) -> float:
    """
    Calculate the mass of the neutral precursor
    Args:
        parsed_pep (list or numba.typed.List of str): the list of amino acids and modified amono acids.
        mass_dict (numba.typed.Dict): key is the amino acid or the modified amino acid, and the value is the mass.
    Returns:
        float: the peptide neutral mass.
    """
    tmass = mass_dict["H2O"]
    for _ in parsed_pep:
        tmass += mass_dict[_]

    return tmass


@njit
def get_fragmass(
    parsed_pep: list,
    mass_dict: numba.typed.Dict
) -> tuple:
    """
    Calculate the masses of the fragment ions
    Args:
        parsed_pep (numba.typed.List of str): the list of amino acids and modified amono acids.
        mass_dict (numba.typed.Dict): key is the amino acid or the modified amino acid, and the value is the mass.
    Returns:
        Tuple[np.ndarray(np.float64), np.ndarray(np.int8)]: the fragment masses and the fragment types (represented as np.int8).
        For a fragment type, positive value means the b-ion, the value indicates the position (b1, b2, b3...); the negative value means
        the y-ion, the absolute value indicates the position (y1, y2, ...).
    """
    n_frags = (len(parsed_pep) - 1) * 2

    frag_masses = np.zeros(n_frags, dtype=np.float64)
    frag_type = np.zeros(n_frags, dtype=np.int8)

    n_frag = 0

    frag_m = mass_dict["Proton"]
    for idx, _ in enumerate(parsed_pep[:-1]):
        frag_m += mass_dict[_]
        frag_masses[n_frag] = frag_m
        frag_type[n_frag] = (idx+1)
        n_frag += 1

    frag_m = mass_dict["Proton"] + mass_dict["H2O"]
    for idx, _ in enumerate(parsed_pep[::-1][:-1]):
        frag_m += mass_dict[_]
        frag_masses[n_frag] = frag_m
        frag_type[n_frag] = -(idx+1)
        n_frag += 1

    return frag_masses, frag_type


def get_frag_dict(
    parsed_pep: list,
    mass_dict: dict
) -> dict:
    """
    Calculate the masses of the fragment ions
    Args:
        parsed_pep (list or numba.typed.List of str): the list of amino acids and modified amono acids.
        mass_dict (numba.typed.Dict): key is the amino acid or the modified amino acid, and the value is the mass.
    Returns:
        dict{str:float}: key is the fragment type (b1, b2, ..., y1, y2, ...), value is fragment mass.
    """
    frag_dict = {}
    frag_masses, frag_type = get_fragmass(parsed_pep, mass_dict)

    for idx, _ in enumerate(frag_masses):

        cnt = frag_type[idx]
        if cnt > 0:
            identifier = 'b'
        else:
            identifier = 'y'
            cnt = -cnt
        frag_dict[identifier+str(cnt)] = _
    return frag_dict


@njit
def calculate_mass(
    mono_mz: float,
    charge: int
) -> float:
    """Calculate the precursor mass from mono mz and charge.
    Args:
        mono_mz (float): mono m/z.
        charge (int): charge.
    Returns:
        float: precursor mass.
    """
    M_PROTON = 1.00727646687
    prec_mass = mono_mz * abs(charge) - charge * M_PROTON
    return prec_mass


@njit
def calculate_mz(
    prec_mass: float,
    charge: int
) -> float:
    """Calculate the precursor mono mz from mass and charge.
    Args:
        prec_mass (float): precursor mass.
        charge (int): charge.
    Returns:
        float: mono m/z.
    """
    M_PROTON = 1.00727646687
    mono_mz = prec_mass / abs(charge) + M_PROTON

    return mono_mz
