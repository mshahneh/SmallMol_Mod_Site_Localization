"""
GNPS Utils - Molecule Utils
---------------------------
This file contains utility functions around ms spectrums

Author: Shahneh
"""

import numpy as np
# from modifinder.utilities.network import get_data
# from modifinder.convert import parse_data_to_universal
from modifinder.classes.Spectrum import Spectrum
import copy

def normalize_peaks(peaks: Spectrum) -> Spectrum:
    """
    l2 normalizes the peaks.
    """
    # l2 normalize the peaks over the intensity
    l2_norm = np.linalg.norm(peaks.intensity)
    new_intensity = [intensity / l2_norm for intensity in peaks.intensity]
    normalized_peaks = Spectrum(peaks)
    normalized_peaks.intensity = new_intensity
    return normalized_peaks


def filter_peaks(peaks, method, variable):
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

    tempPeaks = copy.deepcopy(peaks)
    # call the appropriate function
    if method == "intensity":
        return intensity(tempPeaks, variable)
    elif method == "top_k":
        return top_k(tempPeaks, variable)
    elif method == "both":
        tempPeaks = intensity(tempPeaks, variable)
        return top_k(tempPeaks, 1/variable)
    else:
        return tempPeaks
    

# def get_spectrum(data=None, needs_parse = True, **kwargs):
#     """
#     returns a spectrum from the data

#     Parameters:
#         :data: passed data, can be a USI (str), a dictionary with the necassary keys, a Spectrum object (will return the same), or None
#         :kwargs: keyword arguments, can contain the keys: precursor_mz, precursor_charge, mz, intensity, etc to construct a Spectrum object
#     """

#     if data is None:
#         if needs_parse:
#             data = parse_data_to_universal(kwargs)
#             # import json
#             # print("in get_spectrum\n", json.dumps(data, indent=4))
#         # print("before returning")
#         spec = Spectrum(incoming_data=data)
#         # print("SPEC IS ", spec)
#         return spec

#     elif isinstance(data, str):
#         data = get_data(data)
#         return Spectrum(**data)
    
#     elif isinstance(data, dict):
#         if needs_parse:
#             data = parse_data_to_universal(data)
#         return Spectrum(**data)
    
#     elif isinstance(data, Spectrum):
#         if needs_parse:
#             parsed_data = parse_data_to_universal(data.__dict__)
#             for key in parsed_data:
#                 setattr(data, key, parsed_data[key])
#         return data
    
#     else:
#         raise ValueError("Data type not supported")

