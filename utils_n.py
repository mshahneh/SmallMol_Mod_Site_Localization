import copy
import re
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import collections
import math

SpectrumTuple = collections.namedtuple(
    "SpectrumTuple", ["precursor_mz", "precursor_charge", "mz", "intensity"]
)

def parse_adduct(adduct):
    # adducts = ["+H", "+NH4", "+Na", "+K", "-OH", "-H", "+Cl"]
    adducts = ["+H"]
    acceptedAdducts = ["M" + a for a in adducts]
    if "[" in adduct:
        adduct = adduct.split("[")[1]
        adduct = adduct.split("]")[0]
    if adduct not in acceptedAdducts:
        raise ValueError("Adduct not supported: " + adduct)
    return adduct


def filter_peaks(peaks, method, variable, precursor_mz = None, charge = None):
    """
    Filters the peaks based on the method and variable.
    """
    eps = 0.0001
    def top_k(peaks, k):
        """Filters the peaks by top k peaks."""
        k = int(k)
        filtered_peaks = []
        peaks.sort(key=lambda x: x[1], reverse=True)
        for i in range(min(k, len(peaks))):
            filtered_peaks.append(peaks[i])
        filtered_peaks.sort(key=lambda x: x[0])
        peaks = filtered_peaks
        return filtered_peaks
    
    def intensity(peaks, intensity_ratio_threshold):
        """Filters the peaks by intensity ratio to the maximum peak."""
        filtered_peaks = []
        max_intensity = max([peak[1] for peak in peaks])
        for peak in peaks:
            if peak[1] / max_intensity > intensity_ratio_threshold:
                filtered_peaks.append(peak)
        filtered_peaks.sort(key=lambda x: x[0])
        peaks = filtered_peaks
        return filtered_peaks
    
    # # delete peaks that are bigger than the precursor + charge
    # if precursor_mz != None and charge != None:
    #     peaks = [peak for peak in peaks if peak[0] <= precursor_mz + charge + eps]

    tempPeaks = normalize_peaks(peaks)

    # call the appropriate function
    if method == "intensity":
        return normalize_peaks(intensity(tempPeaks, variable))
    elif method == "top_k":
        return normalize_peaks(top_k(tempPeaks, variable))
    elif method == "both":
        tempPeaks = intensity(tempPeaks, variable)
        return normalize_peaks(top_k(tempPeaks, 1/variable))
    else:
        return tempPeaks

def normalize_peaks(peaks):
    """
    l2 normalizes the peaks.
    """
    # l2 normalize the peaks over the intensity
    l2_norm = np.linalg.norm([peak[1] for peak in peaks])
    normalized_peaks = [(peak[0], peak[1] / l2_norm) for peak in peaks]
    return normalized_peaks


def convert_to_SpectrumTuple(peaks, precursor_mz, precursor_charge):
    """
    Converts the peaks to SpectrumTuple.
    """
    peaks = normalize_peaks(peaks)
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

def find_mz_in_sirius(fragments, search_mz, mz_threshold, ppm_threshold):
    # Binary search in sirius results to find matching mass
    left = 0
    right = len(fragments)

    while left < right:
        mid = (left + right) // 2
        mz = fragments[mid]['mz']
        diff = abs(mz - search_mz)
        if right - left == 1:
            if diff <= mz_threshold and diff <= (ppm_threshold * search_mz / 1e6):
                return left
            else:
                return -1
        elif mz >= search_mz:
            left = mid
        else:
            right = mid

    return -1  # Return -1 if fragment not found


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

    
    modificationSites = []
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in intersect:
            for tempAtom in intersect:
                if mol.GetBondBetweenAtoms(atom.GetIdx(), tempAtom) is not None:
                    modificationSites.append(tempAtom)

    if inParent:
        return modificationSites
    else:
        res = []
        if type(matches[0]) is tuple:
            # multiple matches (substructure is happening multiple times in the parent molecule)
            for match in matches:
                subMatches = list(match)
                temp_res = []
                for atom in modificationSites:
                    idx = subMatches.index(atom) # based on rdkit library, the value in matches array are sorted by index
                    temp_res.append(idx)
                res.append(temp_res)
        else:
            for atom in modificationSites:
                subMatches = list(matches)
                idx = subMatches.index(atom) # based on rdkit library, the value in matches array are sorted by index
                res.append(idx)
        return res
    
def get_modification_graph(main_struct, sub_struct):
    """
        Calculates the substructure difference between main_struct and sub_struct.
        Input:
            main_struct: main molecule
            sub_struct: substructure molecule
        Output:
            count: the graph of the substructure difference, index of the atom of difference that is connected to substructure
    """
    # print(main_struct.GetNumAtoms(), sub_struct.GetNumAtoms())
    if not main_struct.HasSubstructMatch(sub_struct):
        raise ValueError("The substructure is not a substructure of the main structure.")
    atoms_of_modification = []
    atoms_of_substructure = main_struct.GetSubstructMatch(sub_struct)
    for atom in main_struct.GetAtoms():
        if atom.GetIdx() not in atoms_of_substructure:
            atoms_of_modification.append(atom.GetIdx())

    emol = Chem.EditableMol(main_struct)
    for atom in reversed(range(main_struct.GetNumAtoms())):
        if atom not in atoms_of_modification:
            emol.RemoveAtom(atom)
    frag = emol.GetMol()
    if not main_struct.HasSubstructMatch(frag):
        raise ValueError("The substructure is not a substructure of the main structure.")
    modification_atom = calculateModificationSites(main_struct, frag, inParent = False)
    modification_atom = modification_atom[0]

    matches = main_struct.GetSubstructMatch(frag)
    bond_type = Chem.BondType.SINGLE
    for bond in main_struct.GetBonds():
        if bond.GetBeginAtomIdx() in matches and bond.GetEndAtomIdx() not in matches:
            bond_type = bond.GetBondType()
            break
        elif bond.GetEndAtomIdx() in matches and bond.GetBeginAtomIdx() not in matches:
            bond_type = bond.GetBondType()
            break

    return frag, modification_atom, bond_type

def attach_struct(sub_struct, frag, atom, new_modification_atom, bond_type):
    temp_struct = Chem.CombineMols(sub_struct, frag)
    temp_struct = Chem.RWMol(temp_struct)
    temp_struct.AddBond(atom.GetIdx(), new_modification_atom, bond_type)
    # temp_struct.CommitBatchEdit()
    try:
        temp_struct = temp_struct.GetMol()
        # check if mol is valid
        Chem.SanitizeMol(temp_struct)
        modification_atom = temp_struct.GetAtomWithIdx(new_modification_atom)
        return temp_struct
    except:
        temp_struct = Chem.CombineMols(sub_struct, frag)
        temp_struct = Chem.RWMol(temp_struct)
        for bond_type in [Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]:
            temp_struct.AddBond(atom.GetIdx(), new_modification_atom, bond_type)
            temp_struct.CommitBatchEdit()
            try:
                temp_struct = temp_struct.GetMol()
                Chem.SanitizeMol(temp_struct)
                return temp_struct
            except:
                temp_struct = Chem.RWMol(temp_struct)
                temp_struct.RemoveBond(atom.GetIdx(), new_modification_atom)
                temp_struct.CommitBatchEdit()
    
    return None

def generate_possible_stuctures(main_struct, sub_struct):
    def fix_valence(atom, sub_struct):
        if atom.GetNumExplicitHs() > 0:
            temp_struct = Chem.RWMol(sub_struct)
            # remove one explicit H
            neighbors = atom.GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetSymbol() == "H":
                    temp_struct.RemoveBond(atom.GetIdx(), neighbor.GetIdx())
                    temp_struct.CommitBatchEdit()
                    try:
                            temp_struct = temp_struct.GetMol()
                            Chem.SanitizeMol(temp_struct)
                            return temp_struct
                    except:
                        pass
        
        else:
            # if atom has double bond, try to remove one double bond
            neighbors = atom.GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetSymbol() == "H":
                    continue
                temp_struct = Chem.RWMol(sub_struct)
                # remove one double bond
                for bond in temp_struct.GetBonds():
                    if bond.GetBeginAtomIdx() == atom.GetIdx() and bond.GetEndAtomIdx() == neighbor.GetIdx() and bond.GetBondType() == Chem.BondType.DOUBLE:
                        temp_struct.RemoveBond(atom.GetIdx(), neighbor.GetIdx())
                        temp_struct.AddBond(atom.GetIdx(), neighbor.GetIdx(), Chem.BondType.SINGLE)
                        temp_struct.CommitBatchEdit()
                        try:
                            temp_struct = temp_struct.GetMol()
                            Chem.SanitizeMol(temp_struct)
                            return temp_struct
                        except:
                            pass
        
        return None
    
    def _try(sub_struct, frag, atom, new_modification_atom, bond_type):
        temp_struct = attach_struct(sub_struct, frag, atom, new_modification_atom, bond_type)
        if temp_struct is not None:
            # check if the weight of temp_struct is the same as main_struct
            if abs(Descriptors.ExactMolWt(temp_struct) - Descriptors.ExactMolWt(main_struct)) < 0.1:
                return temp_struct
                structs.append((atom.GetIdx(),temp_struct))
        else:
            temp_struct = fix_valence(atom, sub_struct)
            if temp_struct is not None:
                temp_struct = attach_struct(temp_struct, frag, atom, new_modification_atom, bond_type)
                if temp_struct is not None:
                    if abs(Descriptors.ExactMolWt(temp_struct) - Descriptors.ExactMolWt(main_struct)) < 0.1:
                        return temp_struct
                        
                
    frag, modification_atom, bond_type = get_modification_graph(main_struct, sub_struct)
    structs = []
    new_modification_atom = modification_atom + sub_struct.GetNumAtoms()
    # print(new_modification_atom, modification_atom, sub_struct.GetNumAtoms())
    for atom in sub_struct.GetAtoms():
        if atom.GetSymbol() == "H":
            continue
        temp_struct = _try(sub_struct, frag, atom, new_modification_atom, bond_type)
        if temp_struct is not None:
            structs.append((atom.GetIdx(),temp_struct))
        else:
            for additional_type in [Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]:
                temp_struct = _try(sub_struct, frag, atom, new_modification_atom, additional_type)
                if temp_struct is not None:
                    structs.append((atom.GetIdx(),temp_struct))
                    break
                
    return structs


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

def is_submolecule(sub_formula, target_formula, ignore_H = False):
    # Parse the atom counts of the sub-molecule and target molecule
    sub_atom_counts = parse_molecular_formula(sub_formula)
    target_atom_counts = parse_molecular_formula(target_formula)

    # Check if every atom in sub-molecule is in target molecule and has less or equal count
    for element, count in sub_atom_counts.items():
        if element == 'H' and ignore_H:
            continue
        if element not in target_atom_counts or target_atom_counts[element] < count:
            return False

    return True

def add_adduct_to_formula(formula, adduct):
    if adduct != "M+H":
        raise ValueError("Only H adduct is supported.")
    else:
        adduct = copy.deepcopy(adduct)
        adduct = "H"
    # add one H to the formula
    pattern = r"([A-Z][a-z]*)(\d*)"
    matches = re.findall(pattern, formula)  # Find all matches in the formula
    formula = ""
    Found = False
    for match in matches:
        if match[0] == adduct:
            count = match[1]
            if count:
                count = int(count)
            else:
                count = 1
            count += 1
            formula += adduct + str(count)
            Found = True
        else:
            formula += match[0]
            if match[1]:
                formula += match[1]
    if not Found:
        formula += adduct
    
    return formula

def remove_adduct_from_formula(formula, adduct):
    if adduct != "M+H":
        raise ValueError("Only H adduct is supported.")
    else:
        adduct = copy.deepcopy(adduct)
        adduct = "H"
    
    # remove one H from the formula
    pattern = r'([A-Z][a-z]*)(\d*)'
    matches = re.findall(pattern, formula)  # Find all matches in the formula
    formula = ""
    for match in matches:
        if match[0] == "H":
            count = match[1]
            if count:
                count = int(count)
            else:
                count = 1
            count -= 1
            if count > 1:
                formula += "H" + str(count)
            elif count == 1:
                formula += "H"
        else:
            formula += match[0]
            if match[1]:
                formula += match[1]
    
    return formula

def entropy(probabilities):
    # if probabilities is not numpy array, convert it to numpy array
    probabilities = np.array(probabilities)
    # if probabilities is not float, convert it to float
    probabilities = probabilities.astype(float)
    
    if len(probabilities) == 0:
        return 1
    if min(probabilities) < 0:
        probabilities = probabilities - min(probabilities)
    regulator = 1e-8
    probabilities = probabilities + regulator
    probabilities = probabilities / np.sum(probabilities)
    H_max = np.log(len(probabilities))
    H = abs(np.sum(probabilities * np.log(probabilities)))
    # print(H, probabilities, H_max, H/H_max)
    return H/H_max
