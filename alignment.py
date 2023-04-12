import numpy as np
from typing import List, Tuple
import collections
import utils as utils
import handle_network as hn
import json

SpectrumTuple = collections.namedtuple(
    "SpectrumTuple", ["precursor_mz", "precursor_charge", "mz", "intensity"]
)


def norm_intensity(intensity):
    return np.copy(intensity)/np.linalg.norm(intensity)


def _cosine_fast(
    spec: SpectrumTuple,
    spec_other: SpectrumTuple,
    fragment_mz_tolerance: float,
    allow_shift: bool,
) -> Tuple[float, List[Tuple[int, int]]]:
    precursor_charge = max(spec.precursor_charge, 1)
    precursor_mass_diff = (spec.precursor_mz - spec_other.precursor_mz) * precursor_charge
    # Only take peak shifts into account if the mass difference is relevant.
    num_shifts = 1
    if allow_shift and abs(precursor_mass_diff) >= fragment_mz_tolerance:
        num_shifts += precursor_charge
    other_peak_index = np.zeros(num_shifts, np.uint16)
    mass_diff = np.zeros(num_shifts, np.float32)
    for charge in range(1, num_shifts):
        mass_diff[charge] = precursor_mass_diff / charge

    # Find the matching peaks between both spectra.
    peak_match_scores, peak_match_idx = [], []
    for peak_index, (peak_mz, peak_intensity) in enumerate(
        zip(spec.mz, spec.intensity)
    ):
        # Advance while there is an excessive mass difference.
        for cpi in range(num_shifts):
            while other_peak_index[cpi] < len(spec_other.mz) - 1 and (
                peak_mz - fragment_mz_tolerance
                > spec_other.mz[other_peak_index[cpi]] + mass_diff[cpi]
            ):
                other_peak_index[cpi] += 1
                
        # Match the peaks within the fragment mass window if possible.
        for cpi in range(num_shifts):
            index = 0
            other_peak_i = other_peak_index[cpi] + index
            while (
                other_peak_i < len(spec_other.mz)
                and abs(
                    peak_mz - (spec_other.mz[other_peak_i] + mass_diff[cpi])
                )
                <= fragment_mz_tolerance
            ):
                peak_match_scores.append(
                    peak_intensity * spec_other.intensity[other_peak_i]
                )
                peak_match_idx.append((peak_index, other_peak_i))
                index += 1
                other_peak_i = other_peak_index[cpi] + index

    score, peak_matches = 0.0, []
    if len(peak_match_scores) > 0:
        # Use the most prominent peak matches to compute the score (sort in
        # descending order).
        peak_match_scores_arr = np.asarray(peak_match_scores)
        peak_match_order = np.argsort(peak_match_scores_arr)[::-1]
        peak_match_scores_arr = peak_match_scores_arr[peak_match_order]
        peak_match_idx_arr = np.asarray(peak_match_idx)[peak_match_order]
        peaks_used, other_peaks_used = set(), set()
        for peak_match_score, peak_i, other_peak_i in zip(
            peak_match_scores_arr,
            peak_match_idx_arr[:, 0],
            peak_match_idx_arr[:, 1],
        ):
            if (
                peak_i not in peaks_used
                and other_peak_i not in other_peaks_used
            ):
                score += peak_match_score
                # Save the matched peaks.
                peak_matches.append((peak_i, other_peak_i))
                # Make sure these peaks are not used anymore.
                peaks_used.add(peak_i)
                other_peaks_used.add(other_peak_i)

    return score, peak_matches

def handle_usi(molUsi):
    """
    Handles the case where molData and modifData are USIs
    """
    molData = hn.getDataFromUsi(molUsi)
    return molData

def handle_spectrumTuple(molData):
    """
    Handles the case where molData and modifData are SpectrumTuples
    """
    molData = utils.SpectrumTuple_to_dict(molData)
    return molData

def handle_dict(molData):
    """
    Handles the case where molData and modifData are dicts
    """
    molData['peaks'] = json.loads(molData['peaks_json'])
    molData['precursor_mz'] = float(molData['Precursor_MZ'])
    molData['precursor_charge'] = int(molData['Charge'])
    return molData


def handle_alignment(molData, modifData, args={'filter_peaks_method':"top_k", 'filter_peaks_variable':50, 'mz_tolerance':0.05}):

    if args is None:
        raise ValueError("args must be provided")

    if type(molData) == str:
        molData = handle_usi(molData)
    elif type(molData) == SpectrumTuple:
        molData = handle_spectrumTuple(molData)
    elif type(molData) == dict:
        molData = handle_dict(molData)
    else:
        raise TypeError("molData must be of type str, SpectrumTuple, or dict")
    
    if type(modifData) == str:
        modifData = handle_usi(modifData)
    elif type(modifData) == SpectrumTuple:
        modifData = handle_spectrumTuple(modifData)
    elif type(modifData) == dict:
        modifData = handle_dict(modifData)
    else:
        raise TypeError("modifData must be of type str, SpectrumTuple, or dict")
    
    print (molData['peaks'])
    temp1 = utils.filter_peaks(molData['peaks'], args['filter_peaks_method'], args['filter_peaks_variable'])
    molData['peaks'] = temp1
    modifData['peaks'] = utils.filter_peaks(modifData['peaks'], args['filter_peaks_method'], args['filter_peaks_variable'])
    cosine, matchedPeaks = _cosine_fast(utils.convert_to_SpectrumTuple(molData['peaks'], molData["precursor_mz"], molData["precursor_charge"]), 
                                utils.convert_to_SpectrumTuple(modifData['peaks'], modifData["precursor_mz"], modifData["precursor_charge"]),
                                args['mz_tolerance'], allow_shift=True)
    
    res =  {'molData': molData, 'modifData': modifData, 'cosine': cosine, 'matchedPeaks': matchedPeaks}
    return res

    

