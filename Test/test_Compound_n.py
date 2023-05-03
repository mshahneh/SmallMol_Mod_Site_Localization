import sys
sys.path.append('..')
import unittest
from Compound_n import Compound
import json
import os
import pickle
import copy

class TestMoleculeClass(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestMoleculeClass, self).__init__(*args, **kwargs)
        self.library = "BERKELEY-LAB"
        with open(os.path.join("../data/libraries",self.library,"data_dict_filtered.pkl"), "rb") as f:
            self.data_dict_filtered = pickle.load(f)
        self.helperDirectory = os.path.join("../data/libraries",self.library,"nf_output/fragmentationtrees/")
        self.args = {"ppm": 1.01, "mz_tolerance": 0.1}
    
    # TODO: test if can rebuild compund from repr
    def test_repr(self):
        pass
        # compound = Compound(self.compound)
        # stringified = json.loads(str(compound))
        # for key in self.compound:
        #     self.assertTrue(key in stringified)
        #     self.assertEqual(stringified[key], self.compound[key])
    
    def test_generate_peak_to_fragment_map(self):
        cound_not_found = 0
        test_count = 100
        compound_list = list(self.data_dict_filtered.keys())
        for i, compound in enumerate(compound_list):
            if (i > test_count):
                break
            try:
                c = Compound(self.data_dict_filtered[compound], args=self.args)
            except ValueError as e:
                # check if the error is due to an unsupported adduct
                if (str(e).startswith("Adduct not supported:")):
                    continue
                else:
                    raise e
            if (c.Adduct != "M+H"):
                continue
            count = 0
            for i, peak in enumerate(c.peaks):
                if (len (c.peak_fragments_map[i]) == 0):
                    count += 1
            if (len(c.peaks) > 0):
                self.assertFalse(count == len(c.peaks))
            
        # TODO: test if the peak_fragments_map is correct and complete
    
    def test_apply_sirius(self):
        # print count of files in helperDirectory
        file_not_found = 0
        test_count = 100
        test_done = 0
        for compound in self.data_dict_filtered:
            if test_done == test_count:
                break
            try:
                c = Compound(self.data_dict_filtered[compound], args=self.args)
            except ValueError as e:
                # check if the error is due to an unsupported adduct
                if (str(e).startswith("Adduct not supported:")):
                    continue
                else:
                    raise e
            if (c.Adduct != "M+H"):
                continue
            
            before = copy.deepcopy(c.peak_fragments_map)
            try:
                with open(os.path.join(self.helperDirectory, compound + "_fragmentationtree.json")) as f:
                    molSirius = json.load(f)
            except FileNotFoundError:
                file_not_found += 1
                continue
            c.apply_sirius(molSirius)
            for i, peak in enumerate(c.peaks):
                self.assertTrue(len(c.peak_fragments_map[i]) <= len(before[i]))
                if (len(c.peak_fragments_map[i]) == 0):
                    self.assertTrue(len(before[i]) == 0)
            
            test_done += 1
        if (file_not_found > test_count/10):
            print("WARNING: " + str(file_not_found) + " files were not found in " + self.helperDirectory)

if __name__ == '__main__':
    unittest.main()