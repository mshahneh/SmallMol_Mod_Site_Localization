import re
import numpy as np
from rdkit import Chem

def calculateModificationSites(mol, substructure, inParent = True):
    """
        Calculates the number of modification sites to get mol from substructure.
        Input:
            mol1: first molecule
            substructure: substructure molecule
            inParent: bool, if True, the modification sites are given in the parent molecule, if False, the modification sites are given in the substructure
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

def peaks_string_to_list(peaks_string):
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
    filtered_peaks.sort(key=lambda x: x[0])
    peaks = filtered_peaks
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
    peaks = filtered_peaks
    return filtered_peaks

def filter_peaks(peaks, method, variable):
    """
    Filters the peaks based on the method and variable.
    """
    import copy
    tempPeaks = copy.deepcopy(peaks)
    if method == "intensity":
        return filter_peaks_by_intensity(tempPeaks, variable)
    elif method == "top_k":
        return filter_peaks_by_top_k(tempPeaks, variable)
    else:
        return tempPeaks


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

def SpectrumTuple_to_dict(spectrum_tuple):
    """
    Converts the SpectrumTuple to a dictionary.
    """
    res = {}
    res['precursor_charge'] = spectrum_tuple.precursor_charge
    res['precursor_mz'] = spectrum_tuple.precursor_mz
    res['peaks'] = []
    for i in range(len(spectrum_tuple.mz)):
        res['peaks'].append((spectrum_tuple.mz[i], spectrum_tuple.intensity[i]))
    return res

def normalize_peaks(peaks):
    """
    l2 normalizes the peaks.
    """
    # l2 normalize the peaks over the intensity
    l2_norm = np.linalg.norm([peak[1] for peak in peaks])
    peaks = [(peak[0], peak[1] / l2_norm) for peak in peaks]
    return peaks

def is_substruct_substruct(mol1, mol2, indx1, indx2):
    """
    Checks if the substructure mol1[indx1] is a substructure of the substructure mol2[indx2].
    """
    emol = Chem.EditableMol(mol1)
    for atom in reversed(range(mol1.GetNumAtoms())):
        if atom not in indx1:
            emol.RemoveAtom(atom)
    frag1 = emol.GetMol()
    # make all bonds single
    for bond in frag1.GetBonds():
        bond.SetBondType(Chem.rdchem.BondType.SINGLE)

    emol = Chem.EditableMol(mol2)
    for atom in reversed(range(mol2.GetNumAtoms())):
        if atom not in indx2:
            emol.RemoveAtom(atom)
    frag2 = emol.GetMol()
    # make all bonds single
    for bond in frag2.GetBonds():
        bond.SetBondType(Chem.rdchem.BondType.SINGLE)

    return frag2.HasSubstructMatch(frag1)


def parse_molecular_formula(formula):
    # Define the regular expression pattern to match element symbols and their counts
    pattern = r'([A-Z][a-z]*)(\d*)'
    matches = re.findall(pattern, formula)  # Find all matches in the formula

    # Create a dictionary to store element symbol and count pairs
    atom_counts = {}
    for match in matches:
        element = match[0]
        element = element.capitalize()
        count = match[1]
        if count:
            count = int(count)
        else:
            count = 1
        atom_counts[element] = count

    return atom_counts

def is_submolecule(sub_formula, target_formula):
    # Parse the atom counts of the sub-molecule and target molecule
    sub_atom_counts = parse_molecular_formula(sub_formula)
    target_atom_counts = parse_molecular_formula(target_formula)

    # Check if every atom in sub-molecule is in target molecule and has less or equal count
    for element, count in sub_atom_counts.items():
        if element == 'H':
            continue
        if element not in target_atom_counts or target_atom_counts[element] < count:
            return False

    return True