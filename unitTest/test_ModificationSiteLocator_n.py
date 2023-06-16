import sys
sys.path.append('..')
import unittest
from Compound_n import Compound
from ModificationSiteLocator import ModificationSiteLocator
import json
import os
import pickle
import requests
import handle_network as hn


class TestModificationSiteLocator(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestModificationSiteLocator, self).__init__(*args, **kwargs)
        with open(os.path.join("../data/libraries","GNPS-MSMLS","data_dict_filtered.pkl"), "rb") as f:
            self.data_dict_filtered = pickle.load(f)
        with open(os.path.join("../data/libraries", "GNPS-MSMLS", "matches.pkl"), "rb") as f:
            self.matches = pickle.load(f)
    
    # TODO: test the __repr__ function
    def test_repr(self):
        pass
    
    def test_alignment(self):
        test_count = 20
        for match in self.matches[1]:
            try:
                c2, c1 = match
                main_compound = Compound(self.data_dict_filtered[c1], args={"filter_peaks_method":"none"})
                modified_compound = Compound(self.data_dict_filtered[c2], args={"filter_peaks_method":"none"})
                site_locator = ModificationSiteLocator(main_compound, modified_compound, {"mz_tolerance": 0.1})
            except ValueError:
                continue
            self.assertTrue(site_locator.cosine > 0.0 and site_locator.cosine <= 1.0)
            self.assertTrue(len(site_locator.matched_peaks) > 0)
            try:
                matchedPeaks = hn.getMatchedPeaks(hn.generate_usi(c1, "GNPS-MSMLS"), hn.generate_usi(c2, "GNPS-MSMLS"))
            except:
                continue
            try:
                self.assertTrue(site_locator.cosine > matchedPeaks["cosine"]*0.9)
                self.assertTrue(abs(len(site_locator.matched_peaks) - matchedPeaks["n_peak_matches"])/len(site_locator.matched_peaks) < 0.2)
            except AssertionError:
                print("AssertionError: ", c1, c2)
                print("site_locator.cosine: ", site_locator.cosine)
                print("matchedPeaks['cosine']: ", matchedPeaks["cosine"])
                print("len(site_locator.matched_peaks): ", len(site_locator.matched_peaks), site_locator.matched_peaks)
                print("matchedPeaks['n_peak_matches']: ", len(matchedPeaks['peak_matches']), matchedPeaks['peak_matches'])
                raise AssertionError
            test_count -= 1
            if test_count == 0:
                break
    
    # TODO: implement the following function
    def test_find_contributions(self):
        # should only consider the peaks passed in
        # list size should be the same as the number of atoms
        # should contain all the correct fragments
        pass

    # TODO: implement the following function
    def calculate_contributions(self):
        # check PPO
        # check CO
        pass

    # TODO: implement the following function
    def test_generate_probabilities(self):
        pass

    # TODO: implement the following function
    def test_get_structures_by_peak_id(self):
        pass


if __name__ == '__main__':
    unittest.main()