"""
This module provides utilities for handling GNPS (Global Natural Products Social) data types, 
including functions to convert data to and from SpectrumTuple objects, and to parse data into 
a universal format.
Classes:
    SpectrumTuple: A named tuple representing a spectrum with precursor m/z, precursor charge, 
                   m/z values, and intensity values.
Functions:
    convert_to_SpectrumTuple(peaks: list, precursor_mz: float, precursor_charge: int) -> SpectrumTuple:
    convert_to_SpectrumTuple_seprated(mz: list, intensity: list, precursor_mz: float, precursor_charge: int) -> SpectrumTuple:
        Converts separate lists of m/z and intensity values to a SpectrumTuple.
    convert_to_universal_key(key: str) -> str:
        Converts different types of keys to universal keys.
    parse_data_to_universal(data: dict) -> dict:
        Parses the data to a universal format.
    Convert_SpectrumTuple_to_peaks(spectrum: SpectrumTuple) -> list:
        Converts a SpectrumTuple to a list of peaks.
"""

import collections
import json

adduct_mapping = {'M+H': '[M+H]+',
'[M+H]': '[M+H]+',
'[M+H]+': '[M+H]+',
'M+H]': '[M+H]+',
'M+Na': '[M+Na]+',
'[M+Na]': '[M+Na]+',
'[M+Na]+': '[M+Na]+',
'2M+Na': '[2M+Na]+',
'M2+Na': '[2M+Na]+',
'[2M+Na]+': '[2M+Na]+',
'[2M+Na]': '[2M+Na]+',
'M+K': '[M+K]+',
'[M+K]': '[M+K]+',
'[M+K]+': '[M+K]+',
'[2M+K]+': '[2M+K]+',
'2M+K': '[2M+K]+',
'[2M+K]': '[2M+K]+',
'M+H-H20': '[M-H2O+H]+',
'M+H-H2O': '[M-H2O+H]+',
'[M-H2O+H]+': '[M-H2O+H]+',
'M-H20+H': '[M-H2O+H]+',
'[M+H-H2O]+': '[M-H2O+H]+',
'M-H2O+H': '[M-H2O+H]+',
'M+H-2H2O': '[M-2H2O+H]+',
'M-2H2O+H': '[M-2H2O+H]+',
'[M-2H2O+H]+': '[M-2H2O+H]+',
'M-2(H2O)+H': '[M-2H2O+H]+',
'2M+Na-2H': '[2M-2H+Na]-',
'2M-2H+Na': '[2M-2H+Na]-',
'M-H': '[M-H]-',
'[M-H]': '[M-H]-',
'[M-H]-': '[M-H]-',
'M-H-': '[M-H]-',
'M-H1': '[M-H]-',
'3M+Na': '[3M+Na]+',
'[3M+Na]+': '[3M+Na]+',
'[M]+': '[M]+',
'M+': '[M]+',
'M-e': '[M]+',
'M2+H': '[2M+H]+',
'2M+H': '[2M+H]+',
'[2M+H]+': '[2M+H]+',
'[2M+H]': '[2M+H]+',
'[M+2H]': '[M+2H]2+',
'[M+2H]2+': '[M+2H]2+',
'M+2H]': '[M+2H]2+',
'M+2H+2': '[M+2H]2+',
'M+2H': '[M+2H]2+',
'M+acetate': '[M+CH3COOH-H]-',
'M+CH3COOH-H': '[M+CH3COOH-H]-',
'M+CH3COO': '[M+CH3COOH-H]-',
'M+ACN+H': '[M+CH3CN+H]+',
'[M+ACN+H]+': '[M+CH3CN+H]+',
'[M+H+CH3CN]': '[M+CH3CN+H]+',
'M+2Na': '[M+2Na]2+',
'M+2Na]': '[M+2Na]2+',
'M+HCOO': '[M+HCOOH-H]-',
'[M-H+HCOOH]': '[M+HCOOH-H]-',
'M+FA-H': '[M+HCOOH-H]-',
'M+formate': '[M+HCOOH-H]-',
'[M+H+HCOOH]': '[M+HCOOH-H]-',
'2M+FA-H': '[2M+HCOOH-H]-',
'[2M-H+HCOOH]': '[2M+HCOOH-H]-',
'M+NH4': '[M+NH3+H]+',
'[M+NH4]+': '[M+NH3+H]+',
'[M+NH4]': '[M+NH3+H]+',
'2M+Hac-H': '[2M+CH3COOH-H]-',
'2M-H': '[2M-H]-',
'[2M-H]': '[2M-H]-',
'2M+NH4': '[2M+NH3+H]+',
'[2M+NH4]+': '[2M+NH3+H]+',
'[2M+NH4]': '[2M+NH3+H]+',
'[2M+Ca]2+': '[2M+Ca]2+',
'[M+Ca]2+': '[M+Ca]2+',
'[3M+Ca]2+': '[3M+Ca]2+',
'[2M+Ca-H]+': '[2M-H+Ca]+',
'[2M-H2O+H]+': '[2M-H2O+H]+',
'[4M+Ca]2+': '[4M+Ca]2+',
'[3M+NH4]+': '[3M+NH3+H]+',
'3M+NH4': '[3M+NH3+H]+',
'[2M-2H2O+H]+': '[2M-2H2O+H]+',
'[M+ACN+NH4]+': '[M+CH3CN+NH3+H]+',
'[5M+Ca]2+': '[5M+Ca]2+',
'[3M+K]+': '[3M+K]+',
'[3M+Ca-H]+': '[3M-H+Ca]2+',
'[M-H+2Na]+': '[M-H+2Na]+',
'M-H+2Na': '[M-H+2Na]+',
'[M-3H2O+H]+': '[M-3H2O+H]+',
'M-3H2O+H': '[M-3H2O+H]+',
'[M-3H2O+2H]2+': '[M-3H2O+2H]2+',
'[M-2H2O+2H]2+': '[M-2H2O+2H]2+',
'[M-4H2O+H]+': '[M-4H2O+H]+',
'[M-5H2O+H]+': '[M-5H2O+H]+',
'[M+Ca-H]+': '[M+Ca-H]+',
'[2M-H+2Na]+': '[2M-H+2Na]+',
'[2M-3H2O+H]+': '[2M-3H2O+H]+',
'[M+H+Na]2+': '[M+Na+H]2+',
'[M-2H2O+NH4]+': '[M-2H2O+NH3+H]+',
'[2M-2H+Na]': '[2M-2H+Na]-',
'[M-H+CH3OH]': '[M+CH3OH-H]-',
'M+MeOH-H': '[M+CH3OH-H]-',
'M-H2O-H': '[M-H2O-H]-',
'[M-H-H2O]': '[M-H2O-H]-',
'M+Cl-': '[M+Cl]-',
'M+Cl': '[M+Cl]-',
'[M+Cl]': '[M+Cl]-',
'M+K-2H': '[M-2H+K]-',
'[M-2H+K]': '[M-2H+K]-',
'M-2H]': '[M-2H]2-',
'M-2H': '[M-2H]2-',
'M-2H-': '[M-2H]2-',
'M+Na-2H': '[M-2H+Na]-',
'[M-2H+Na]': '[M-2H+Na]-',
'M+Br': '[M+Br]-',
'3M-H': '[3M-H]-',
'[3M-H]': '[3M-H]-',
'[M+H+CH3OH]': '[M+CH3OH+H]+',
'M+CH3OH+H': '[M+CH3OH+H]+',
'[2M+H+CH3CN]': '[2M+CH3CN+H]+',
'M-CO2-H': '[M-CO2-H]-',
'[2M-2H+K]': '[2M-2H+K]-',
'2M+K-2H': '[2M-2H+K]-',
'[M+Na+CH3CN]': '[M+CH3CN+Na]+',
'M-H2+H': '[M-H2+H]-',
'M-H+Cl]': '[M-H+Cl]2-',
'M-H+Cl': '[M-H+Cl]2-',
'3M+H': '[3M+H]+',
'[3M+H]': '[3M+H]+',
'M+H-NH3': '[M-NH3+H]+',
'M-NH3+H': '[M-NH3+H]+',
'M-H+C2H2O': '[M+C2H2O-H]-',
'M+H-C2H2O': '[M+C2H2O-H]-',
'M-H+CH2O2': '[M+CH2O2-H]-',
'M+CH2O2-H': '[M+CH2O2-H]-',
'M+TFA-H': '[M+C2HF3O2-H]-',
'M-C2HF3O2-H': '[M+C2HF3O2-H]-',
'[M]1+': '[M]1+'}

SpectrumTuple = collections.namedtuple(
    "SpectrumTuple", ["precursor_mz", "precursor_charge", "mz", "intensity"]
)


def convert_to_SpectrumTuple(peaks: list, precursor_mz: float, precursor_charge: int) -> SpectrumTuple:
    """
    Converts a list of peaks to a SpectrumTuple.

    Args:
        :peaks (list): A list of tuples (mz, intensity) or a list of lists [mz, intensity].
        :precursor_mz (float): The precursor m/z value.
        :precursor_charge (int): The precursor charge value.
    Returns:
        :SpectrumTuple: An instance of SpectrumTuple containing the provided peaks and precursor information.
    Raises:
        :ValueError: If peaks is not a list of tuples or lists.
    """
    if not peaks:
        return None
    if len(peaks) > 0 and (not isinstance(peaks[0], tuple) and not isinstance(peaks[0], list)):
        raise ValueError("Peaks should be a list of tuples (mz, intensity)")
    
    res = {}
    res['precursor_charge'] = precursor_charge
    res['precursor_mz'] = precursor_mz
    res['mz'] = []
    res['intensity'] = []
    for peak in peaks:
        res['mz'].append(peak[0])
        res['intensity'].append(peak[1])
    
    return SpectrumTuple(**res)

def convert_to_SpectrumTuple_seprated(mz: list, intensity: list, precursor_mz: float, precursor_charge: int) -> SpectrumTuple:
    """
    Converts separate lists of m/z and intensity values to a SpectrumTuple.

    Args:
        :mz (list): List of mz values.
        :intensity (list): List of intensity values.
        :precursor_mz (float): Precursor m/z.
        :precursor_charge (int): Precursor charge.
    Returns:
        :SpectrumTuple: A SpectrumTuple object containing the provided mz, intensity, precursor_mz, and precursor_charge.
    Raises:
        :ValueError: If the length of mz and intensity lists are not the same.
    """
    if not mz:
        return None
    if len(mz) != len(intensity):
        raise ValueError("Length of mz and intensity should be same")
    
    res = {}
    res['precursor_charge'] = precursor_charge
    res['precursor_mz'] = precursor_mz
    res['mz'] = mz
    res['intensity'] = intensity
    
    return SpectrumTuple(**res)

def convert_to_universal_key(key: str) -> str:
    """
    Convert different types of keys to universal keys.
    This function standardizes various key names to a universal format. 
    It currently supports conversion for the following keys:
    - "precursor_mz" to "Precursor_MZ"
    - "smiles" or "SMILES" to "Smiles"
    - "charge" to "Charge"
    - "adduct" to "Adduct"

    Args:
        :key (str): The key to be converted.
    
    Returns:
        :str: The converted key.
    """

    if key == "precursor_mz":
        return "Precursor_MZ"
    if key == "smiles" or key == "SMILES":
        return "Smiles"
    if key == "charge":
        return "Charge"
    if key == "adduct":
        return "Adduct"
    # TODO: Add more keys
    return key
    
def parse_data_to_universal(data):
    """
    Parse the data to a universal format.

    This function takes a dictionary of data and converts it into a universal format.
    It processes specific keys like "peaks_json" and "Charge" differently, and attempts
    to convert other values to floats. If the conversion to float is successful and the
    key is "Charge", it further converts the value to an integer.

    Args:
        :data (dict): The input data dictionary to be parsed.

    Returns:
        :dict: A dictionary with keys converted to a universal format and values processed
              accordingly.
    """

    res = {}
    for key, value in data.items():
        converted_key = convert_to_universal_key(key)
        if key == "peaks_json":
            res['peaks'] = json.loads(value)
        elif converted_key == "Adduct":
            res[converted_key] = adduct_mapping.get(value, value)
        else:
            try:
                value = float(value)
                if key == "Charge":
                    value = int(value)
            except:
                pass
            res[converted_key] = value
    return res

def Convert_SpectrumTuple_to_peaks(spectrum: SpectrumTuple):
    """
    Convert a SpectrumTuple to a list of peaks.

    This function takes a SpectrumTuple object, which contains m/z (mass-to-charge ratio) and intensity values,
    and converts it into a list of tuples where each tuple represents a peak with its corresponding m/z and intensity.

    Args:
        :spectrum (SpectrumTuple): The SpectrumTuple object containing m/z and intensity values.

    Returns:
        :list of tuple: A list of tuples where each tuple contains an m/z value and its corresponding intensity.
    """
    peaks = []
    for mz, intensity in zip(spectrum.mz, spectrum.intensity):
        peaks.append((mz, intensity))
    return peaks