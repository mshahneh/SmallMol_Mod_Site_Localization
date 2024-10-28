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
        if key == "peaks_json":
            res['peaks'] = json.loads(value)
        else:
            try:
                value = float(value)
                if key == "Charge":
                    value = int(value)
            except:
                pass
            res[convert_to_universal_key(key)] = value
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