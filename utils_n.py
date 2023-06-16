import copy
import re
import numpy as np

def parse_adduct(adduct):
    adducts = ["+H", "+NH4", "+Na", "+K", "-OH", "-H", "+Cl"]
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

def find_mz_in_sirius(fragments, search_mz, threshold):
    # Binary search in sirius results to find matching mass
    left = 0
    right = len(fragments)

    while left < right:
        mid = (left + right) // 2
        mz = fragments[mid]['mz']
        diff = abs(mz - search_mz)
        if right - left == 1:
            if diff <= threshold:
                return left
            else:
                return -1
        elif mz >= search_mz:
            left = mid
        else:
            right = mid

    return -1  # Return -1 if fragment not found

def is_submolecule(sub_formula, target_formula):
    # Parse the atom counts of the sub-molecule and target molecule

    def parse_molecular_formula(formula):
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

    sub_atom_counts = parse_molecular_formula(sub_formula)
    target_atom_counts = parse_molecular_formula(target_formula)

    # Check if every atom in sub-molecule is in target molecule and has less or equal count
    for element, count in sub_atom_counts.items():
        if element == 'H':
            continue
        if element not in target_atom_counts or target_atom_counts[element] < count:
            return False

    return True
