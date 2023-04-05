import os
import json
import requests
import numpy as np

def calculateModificationSites(mol, substructure, inParent = True):
    """
        Calculates the number of modification sites to get mol from substructure.
        Input:
            mol1: first molecule
            substructure: substructure molecule
        Output:
            count: modification sites
    """

    matches = mol.GetSubstructMatch(substructure)
    intersect = set(matches)

    
    modificationSites = set()
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in intersect:
            for tempAtom in intersect:
                if mol.GetBondBetweenAtoms(atom.GetIdx(), tempAtom) is not None:
                    modificationSites.add(tempAtom)

    if inParent:
        return modificationSites
    else:
        res = set()
        for atom in modificationSites:
            if type(matches[0]) is tuple:
                for match in matches:
                    subMatches = list(match)
                    idx = subMatches.index(atom)
                    res.add(idx)
            else:
                subMatches = list(matches)
                idx = subMatches.index(atom)
                res.add(idx)
        return list(res)

def getHitAtomsAndBonds(mol, substructure):
    """
        Returns the atoms and bonds that are hit by the substructure.
        Input:
            mol: molecule
            substructure: substructure molecule
        Output:
            hitAtoms: list of hit atoms
            hitBonds: list of hit bonds
    """
    
    matches = mol.GetSubstructMatches(substructure)
    hitAtoms = []
    hitBonds = []
    if len(matches) == 0:
        return hitAtoms, hitBonds
    
    if type(matches[0]) is tuple:
        for match in matches:
            tempHitAtoms = set()
            tempHitBonds = set()
            for atom in match:
                tempHitAtoms.add(atom)
            for bond in mol.GetBonds():
                if bond.GetBeginAtomIdx() in tempHitAtoms and bond.GetEndAtomIdx() in tempHitAtoms:
                    tempHitBonds.add(bond.GetIdx())
            if len(tempHitAtoms) > 0:
                hitAtoms.append(list(tempHitAtoms))
                hitBonds.append(list(tempHitBonds))
    else:
        hitAtoms.append([])
        hitBonds.append([])
        for atom in matches:
            hitAtoms[0].append(atom)
        for bond in mol.GetBonds():
            if bond.GetBeginAtomIdx() in hitAtoms[0] and bond.GetEndAtomIdx() in hitAtoms[0]:
                hitBonds[0].append(bond.GetIdx())

    return hitAtoms, hitBonds


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
    r = requests.get('https://metabolomics-usi.ucsd.edu/json/mirror/', params=payload,  timeout=5)
    return json.loads(r.text)

def getDataFromUsi(usi):
    url = 'https://metabolomics-usi.ucsd.edu/json/' + "?usi1=" + usi
    r = requests.get(url)
    return json.loads(r.text)

def separateShifted(matchedPeaks, mol1peaks, mol2peaks, eps = 0.1):
    """
    Separates the shifted and unshifted peaks.
    """
    shifted = []
    unshifted = []
    for peak in matchedPeaks:
        if abs(mol1peaks[peak[0]][0] - mol2peaks[peak[1]][0]) > eps:
            shifted.append(peak)
        else:
            unshifted.append(peak)
    return shifted, unshifted

def concert_explaind_peaks_string_to_list(peaks_string):
    """
    Converts the explained peaks string to a list of tuples.
    """
    peaks_list = []
    for peak in peaks_string.split(";"):
        peak = peak.split(":")
        peaks_list.append((float(peak[0]), peak[1]))
    return peaks_list


def filter_peaks_by_intensity(peaks, intensity_ratio_threshold):
    """
    Filters the peaks by intensity ratio to the maximum peak.
    """
    filtered_peaks = []
    max_intensity = max([peak[1] for peak in peaks])
    for peak in peaks:
        if peak[1] / max_intensity > intensity_ratio_threshold:
            filtered_peaks.append(peak)
    return filtered_peaks

def filter_peaks_by_top_k(peaks, k):
    """
    Filters the peaks by top k peaks.
    """
    k = int(k)
    filtered_peaks = []
    peaks.sort(key=lambda x: x[1], reverse=True)
    for i in range(min(k, len(peaks))):
        filtered_peaks.append(peaks[i])
    filtered_peaks.sort(key=lambda x: x[0])
    return filtered_peaks

def filter_peaks(peaks, method, variable):
    """
    Filters the peaks based on the method and variable.
    """
    if method == "intensity":
        return filter_peaks_by_intensity(peaks, variable)
    elif method == "top_k":
        return filter_peaks_by_top_k(peaks, variable)
    else:
        return peaks
    
def convert_to_SpectrumTuple(peaks, precursor_mz, precursor_charge):
    """
    Converts the peaks to SpectrumTuple.
    """
    peaks = normalize_peaks(peaks)
    from alignment import SpectrumTuple
    res = {}
    res['precursor_charge'] = precursor_charge
    res['precursor_mz'] = precursor_mz
    res['mz'] = []
    res['intensity'] = []
    for peak in peaks:
        res['mz'].append(peak[0])
        res['intensity'].append(peak[1])
    
    return SpectrumTuple(**res)

def normalize_peaks(peaks):
    """
    l2 normalizes the peaks.
    """
    # l2 normalize the peaks over the intensity
    l2_norm = np.linalg.norm([peak[1] for peak in peaks])
    peaks = [(peak[0], peak[1] / l2_norm) for peak in peaks]
    return peaks