"""
GNPS Utils - Network Module

This module provides functionality to connect to GNPS and retrieve data.

Author: Shahneh
"""

import json
import requests
from modifinder.utilities.gnps_types import *

def usi_to_accession(usi: str) -> str:
    """
    Get the accession number from a USI
    param usi: str
    return: str if found, None otherwise
    """

    if "accession" in usi:
        return usi.split(":")[-1]
    else:
        if type(usi) == str:
            return usi
        else:
            return None


def accession_to_usi(accession: str) -> str:
    """
    Get the USI from an accession id
    param accession: str
    return: str
    """
    return "mzspec:GNPS:{}:accession:{}".format("GNPS-LIBRARY", accession)


def get_data(identifier: str) -> dict:
    """
    Get data from GNPS, either from USI or Accession. if the identifier points to a known item in gnps,
      it will return the full data, otherwise it will return partial data (ms2 data)
    param identifier: str - USI or Accession
    return: dict - dictionary of data
    """

    if _is_usi(identifier):
        if _is_known(identifier):
            identifier = usi_to_accession(identifier)
        else:
            data = _get_partial_data(identifier)
            data = parse_data_to_universal(data)
            return data

    link = "https://external.gnps2.org/gnpsspectrum?SpectrumID={}".format(identifier)
    res = requests.get(link)
    parsed = res.json()

    data = dict()
    try:
        data.update(parsed['annotations'][0])
    except KeyError:
        pass
    try:
        data.update(parsed['spectruminfo'])
    except KeyError:
        pass
    try:
        data['comments'] = parsed['comments']
    except KeyError:
        pass

    data = parse_data_to_universal(data)

    return data

def get_matched_peaks(identifier1: str, identifier2: str) -> dict:
    """
    runs the gnps modified cosine matching algorithm and returns the matched peaks
    param identifier1: str - USI or Accession
    param identifier2: str - USI or Accession
    return: dict - dictionary of matched peaks
    """
    if not _is_usi(identifier1):
        identifier1 = accession_to_usi(identifier1)
    
    if not _is_usi(identifier2):
        identifier2 = accession_to_usi(identifier2)

    payload = {
        'usi1': identifier1,
        'usi2': identifier2,
     'mz_min': 'None',
     'mz_max':'None',
     'cosine':'shifted',
     'fragment_mz_tolerance':'0.1',
      'grid': 'True'}
    r = requests.get('https://metabolomics-usi.gnps2.org/json/mirror/', params=payload,  timeout=5)
    return json.loads(r.text)


def _get_partial_data(identifier: str) -> dict:
    """
    Get partial data (ms2 data) from USI
    param identifier: str - USI
    return: dict - dictionary of data with keys: precursor_mz, precursor_charge, mz: list, intensity: list
    """
    url = 'https://metabolomics-usi.gnps2.org/json/' + "?usi1=" + identifier
    r = requests.get(url)
    data = json.loads(r.text)
    # change key names from precaursor_mz to Precursor_MZ
    data['Precursor_MZ'] = data.pop('precursor_mz')
    data['Charge'] = data.pop('precursor_charge')

    return data


def _is_usi(identifier: str) -> bool:
    """
    Check if the identifier is a USI
    param identifier: str
    return: bool
    """
    return "mzspec" in identifier


def _is_known(identifier: str) -> bool:
    """
    Check if the identifier is a known identifier in GNPS
    param identifier: str
    return: bool
    """
    return "accession" in identifier

def get_spectrum(identifier) -> SpectrumTuple:
    """
    Get the spectrum from USI or Accession
    :param identifier: str - USI or Accession
    :return: SpectrumTuple
    """
    if isinstance(identifier, SpectrumTuple):
        return identifier
    
    if isinstance(identifier, list):
        return convert_to_SpectrumTuple(identifier, 0, 0)
    
    if isinstance(identifier, dict):
        return convert_to_SpectrumTuple(identifier['peaks'], identifier['Precursor_MZ'], identifier['Charge'])
    
    data = get_data(identifier)
    return convert_to_SpectrumTuple(data['peaks'], data['Precursor_MZ'], data['Charge'])
    