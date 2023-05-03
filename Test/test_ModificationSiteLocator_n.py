import sys
sys.path.append('..')
import unittest
from Compound_n import Compound
from ModificationSiteLocator import ModificationSiteLocator
import json
import os
import pickle
import requests

def getMatchedPeaks(usi1, usi2):
    payload = {
        'usi1': usi1,
        'usi2': usi2,
     'mz_min': 'None',
     'mz_max':'None',
     'annotate_precision': '2',
     'annotation_rotation':'45',
     'max_intensity': '50',
     'cosine':'shifted',
     'fragment_mz_tolerance':'0.1',
    #  'annotate_peaks': 'value3',
      'grid': 'True'}
    r = requests.get('https://metabolomics-usi.ucsd.edu/json/mirror/', params=payload,  timeout=5)
    return json.loads(r.text)


class TestModificationSiteLocator(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestModificationSiteLocator, self).__init__(*args, **kwargs)
        with open(os.path.join("../data/libraries","GNPS-MSMLS","data_dict_filtered.pkl"), "rb") as f:
            data_dict_filtered = pickle.load(f)
        c1 = "CCMSLIB00005464402"
        c2 = "CCMSLIB00005464389"
        self.main_compound = Compound(data_dict_filtered[c1])
        self.modified_compound = Compound(data_dict_filtered[c2])
        self.usi1 = "mzspec:GNPS:GNPS-LIBRARY:"+c1
        self.usi2 = "mzspec:GNPS:GNPS-LIBRARY:"+c2
        self.data_dict_filtered = data_dict_filtered
    
    # TODO: test the __repr__ function
    def test_repr(self):
        pass
    
    # TODO: test alignment with metablomics-usi data 
    def test_alignment(self):
        site_locator = ModificationSiteLocator(self.main_compound, self.modified_compound)
        self.assertTrue(site_locator.cosine > 0.0 and site_locator.cosine <= 1.0)
        self.assertTrue(len(site_locator.matched_peaks) > 0)
        # matchedPeaks = getMatchedPeaks(self.usi1, self.usi2)
        # self.assertEqual(site_locator.cosine, matchedPeaks["cosine"])
        # self.assertEqual(len(site_locator.matched_peaks), len(matchedPeaks["matchedPeaks"]))


if __name__ == '__main__':
    unittest.main()