import sys
sys.path.append('..')
import unittest
import os
import numpy as np
import pickle
from Compound_n import Compound
from ModificationSiteLocator import ModificationSiteLocator
from SiteLocator import SiteLocator
import utils

# library = "GNPS-LIBRARY"
library = "GNPS-MSMLS"
# library = "GNPS-NIH-NATURALPRODUCTSLIBRARY_ROUND2_POSITIVE"

# load data_dict_filtered
with open(os.path.join("../data/libraries",library,"data_dict_filtered.pkl"), "rb") as f:
    data_dict_filtered = pickle.load(f)

# load matches
with open(os.path.join("../data/libraries",library,"matches.pkl"), "rb") as f:
    matches = pickle.load(f)

# load cachedStructures_filtered
with open(os.path.join("../data/libraries",library,"cachedStructures.pkl"), "rb") as f:
    cachedStructures_filtered = pickle.load(f)


args = {
    'filter_peaks_method':"top_k",
      'filter_peaks_variable':60,
        'mz_tolerance':0.05,
        'mz_tolerance':0.05
}


for setting in range(8):
    peak_presence_only = setting % 2 == 0
    consider_intensity = setting // 2 % 2 == 0
    combine = setting // 4 % 2 == 0
    count = 200
    for match in matches[1]:
        m1, m0 = match
        if data_dict_filtered[m0]["Adduct"] != data_dict_filtered[m1]["Adduct"] or data_dict_filtered[m0]["Adduct"] != "M+H":
            continue
        main_compound = Compound(data_dict_filtered[m0], cachedStructures_filtered[m0], args)
        modified_compound = Compound(data_dict_filtered[m1], cachedStructures_filtered[m1], args)
        site_locator = ModificationSiteLocator(main_compound, modified_compound, args)
        new_scores = site_locator.generate_probabilities(shifted_only=(not combine), PPO=peak_presence_only, CI=consider_intensity)
        modification_site = utils.calculateModificationSites(modified_compound.structure, main_compound.structure, False)[0]
        
        molMol = cachedStructures_filtered[m0]
        modifMol = cachedStructures_filtered[m1]
        molSmiles = data_dict_filtered[m0]['Smiles']
        modifSmiles = data_dict_filtered[m1]['Smiles']
        site = SiteLocator(data_dict_filtered[m0], data_dict_filtered[m1], molSmiles, args)

        scores_unshifted, scores_shifted = site.calculate_score(peak_presence_only, consider_intensity)
        old_scores = site.distance_score(scores_unshifted, scores_shifted, combine)
        old_scores = old_scores - np.min(old_scores)
        if np.sum(old_scores) != 0:
            old_scores = old_scores / np.sum(old_scores)
        else:
            old_scores = np.zeros(len(old_scores))

        # if new_scores not equal to old_scores, throw error
        if not np.allclose(new_scores, old_scores):
            print("new_scores: " + str(new_scores))
            print("old_scores: " + str(old_scores))
            print("new_scores - old_scores: " + str(new_scores - old_scores))
            print("m0", m0)
            print("m1", m1)
            raise ValueError("new_scores not equal to old_scores")
        
        else:
            count -= 1
        
        if count == 0:
            break
    print("setting", setting, "passed")

print("All tests passed!")