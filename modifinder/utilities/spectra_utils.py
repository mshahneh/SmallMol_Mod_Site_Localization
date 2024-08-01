import numpy as np
from modifinder.utilities.network import *
from modifinder.utilities.gnps_types import *
import copy

def normalize_peaks(peaks: SpectrumTuple) -> SpectrumTuple:
    """
    l2 normalizes the peaks.
    """
    # l2 normalize the peaks over the intensity
    l2_norm = np.linalg.norm(peaks.intensity)
    new_intensity = [intensity / l2_norm for intensity in peaks.intensity]
    normalized_peaks = SpectrumTuple(
        precursor_mz=peaks.precursor_mz,
        precursor_charge=peaks.precursor_charge,
        mz=peaks.mz,
        intensity=new_intensity,
    )
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