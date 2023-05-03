import sys
sys.path.append('..')
import unittest
from Compound_n import Compound
from utils_n import parse_adduct, filter_peaks, is_submolecule, find_mz_in_sirius
import json
import os
import pickle

class TestUtils(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestUtils, self).__init__(*args, **kwargs)
        with open(os.path.join("../data/libraries","GNPS-MSMLS","data_dict_filtered.pkl"), "rb") as f:
            self.data_dict_filtered = pickle.load(f)
    
    def test_parse_adduct(self):
        adduct = ["+H", "+NH4", "+Na", "+K", "-OH", "-H", "+Cl"]
        acceptedAdducts = ["M" + a for a in adduct]
        acceptedAdducts.extend(["2M" + a for a in adduct])
        for key in self.data_dict_filtered:
            try:
                parsed_adduct = parse_adduct(self.data_dict_filtered[key]["Adduct"])
            except ValueError:
                continue
            self.assertTrue(parsed_adduct in acceptedAdducts)
    
    # TODO: test filter_peaks
    def test_filter_peaks(self):
        pass

    # TODO: test is_submolecule
    def test_is_submolecule(self):
        pass

    # TODO: test find_mz_in_sirius
    def test_find_mz_in_sirius(self):
        pass

if __name__ == '__main__':
    unittest.main()