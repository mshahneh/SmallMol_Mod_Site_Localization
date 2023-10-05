import os
import json
import requests
import numpy as np
from .utils_n import convert_to_SpectrumTuple


def generate_usi(id, library_membership):
    return "mzspec:GNPS:" + library_membership + ":accession:" + id

def getMatchedPeaks(usi1, usi2):
    payload = {
        'usi1': usi1,
        'usi2': usi2,
     'mz_min': 'None',
     'mz_max':'None',
     'annotate_precision': '2',
     'annotation_rotation':'45',
     'max_intensity': '50',
     'cosine':'shifted',
     'fragment_mz_tolerance':'0.1',
    #  'annotate_peaks': 'value3',
      'grid': 'True'}
    r = requests.get('https://metabolomics-usi.gnps2.org/json/mirror/', params=payload,  timeout=5)
    return json.loads(r.text)

def getDataFromUsi(usi):
    url = 'https://metabolomics-usi.gnps2.org/json/' + "?usi1=" + usi
    r = requests.get(url)
    data = json.loads(r.text)
    # change key names from precaursor_mz to Precursor_MZ
    data['Precursor_MZ'] = data.pop('precursor_mz')
    data['Charge'] = data.pop('precursor_charge')

    # if there are no adducts, add a default one
    if 'Adduct' not in data:
        data['Adduct'] = "M+H"

    return data


def generate_usi(id, library_membership):
    return "mzspec:GNPS:" + library_membership + ":accession:" + id

def getMatchedPeaks(usi1, usi2):
    payload = {
        'usi1': usi1,
        'usi2': usi2,
     'mz_min': 'None',
     'mz_max':'None',
     'annotate_precision': '2',
     'annotation_rotation':'45',
     'max_intensity': '50',
     'cosine':'shifted',
     'fragment_mz_tolerance':'0.1',
    #  'annotate_peaks': 'value3',
      'grid': 'True'}
    r = requests.get('https://metabolomics-usi.gnps2.org/json/mirror/', params=payload,  timeout=5)
    return json.loads(r.text)

def usi_to_SpectrumTuple(usi):
    """
    Converts the usi to SpectrumTuple.
    """
    data = getDataFromUsi(usi)
    return convert_to_SpectrumTuple(data['peaks'], data['precursor_mz'], data['precursor_charge'])