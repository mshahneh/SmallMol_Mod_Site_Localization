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
    #  'annotate_precision': '2',
    #  'annotation_rotation':'45',
    #  'max_intensity': '50',
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
    # https://external.gnps2.org/gnpsspectrum?SpectrumID=CCMSLIB00005464251

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

def generate_link_accession(id):
    return "https://external.gnps2.org/gnpsspectrum?SpectrumID={}".format(id)

def get_library_from_accession(id, get_smiles=False):
    import requests
    res = requests.get(generate_link_accession(id))
    parsed = res.json()
    if get_smiles:
        try:
            return parsed['spectruminfo']['library_membership'], parsed['annotations'][0]['Smiles']
        except:
            return parsed['spectruminfo']['library_membership'], None
    else:
        return parsed['spectruminfo']['library_membership']

def create_usi_from_accession(id, library = None):
    if library is None:
        library = get_library_from_accession(id)
    return "mzspec:GNPS:{}:accession:{}".format(library, id)

def create_link_from_accession(id1, id2):
    base = "https://modsitelocalization.gnps2.org/"
    # base = "http://localhost:5000/"
    library1, smiles1 = get_library_from_accession(id1, True)
    library2, smiles2 = get_library_from_accession(id2, True)
    usi1 = create_usi_from_accession(id1, library1)
    usi2 = create_usi_from_accession(id2, library2)
    if smiles2 is None:
        url = base + "?USI1=" + usi1 + "&USI2=" + usi2 + "&SMILES1=" + smiles1
    else:
        url = (
            base
            + "?USI1="
            + usi1
            + "&USI2="
            + usi2
            + "&SMILES1="
            + smiles1
            + "&SMILES2="
            + smiles2
        )
    return url