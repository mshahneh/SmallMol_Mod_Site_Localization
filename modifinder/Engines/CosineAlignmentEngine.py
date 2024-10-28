import numpy as np
from typing import List, Tuple
import json
from modifinder.Engines.engine_abstracts import AlignmentEngine
from modifinder.utilities.gnps_types import SpectrumTuple


def _cosine_fast(
    spec: SpectrumTuple,
    spec_other: SpectrumTuple,
    fragment_mz_tolerance: float,
    fragment_ppm_tolerance: float,
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
                and abs(peak_mz - (spec_other.mz[other_peak_i] + mass_diff[cpi])) <= fragment_mz_tolerance
            ):
                if abs(peak_mz - (spec_other.mz[other_peak_i] + mass_diff[cpi])) <= (fragment_ppm_tolerance * peak_mz / 1e6):
                    peak_match_scores.append(peak_intensity * spec_other.intensity[other_peak_i])
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

class CosineAlignmentEngine(AlignmentEngine):
    def __init__(self):
        pass

    def align(self, network, **kwargs):
        pass

    def single_align(self, SpectrumTuple1: SpectrumTuple,
                      SpectrumTuple2: SpectrumTuple, 
                      fragment_mz_tolerance: float = 0.02, 
                      fragment_ppm_tolerance: float = 100.0):
        """
        Aligns two spectra using cosine similarity and returns the cosine score and the matched peaks.

        Parameters:
            SpectrumTuple1 (SpectrumTuple): First spectrum
            SpectrumTuple2 (SpectrumTuple): Second spectrum
            fragment_mz_tolerance (float): Fragment mz tolerance
            fragment_ppm_tolerance (float): Fragment ppm tolerance

        Returns:
            Tuple[float, List[Tuple[int, int]]]: alignment score and the list of matched peaks, each value in the list is a tuple of the indices of the matched peaks in the two spectra
        """
        cosine, matched_peaks = _cosine_fast(
            SpectrumTuple1, SpectrumTuple2, fragment_mz_tolerance, fragment_ppm_tolerance, True
        )
        return cosine, matched_peaks