import numpy as np
from modifinder.utilities.network import *
from modifinder.utilities.gnps_types import *


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