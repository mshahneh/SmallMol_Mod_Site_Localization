import sys
sys.path.append('..')
import unittest
from Compound_n import Compound
from fragmentation_py import FragmentEngine
import rdkit.Chem as Chem
import json
import os
import pickle
import copy

class TestMoleculeClass(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestMoleculeClass, self).__init__(*args, **kwargs)
        with open('test_samples.pkl', 'rb') as f:
            self.samples = pickle.load(f)
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
        for i, sample in enumerate(self.samples):
            compound = sample['main']
            try:
                c = Compound(compound, args=self.args)
            except ValueError as e:
                # check if the error is due to an unsupported adduct
                if (str(e).startswith("Adduct not supported:")):
                    continue
                else:
                    raise e
            if (c.Adduct != "M+H"):
                continue
            count = 0
            for j, peak in enumerate(c.peaks):
                if (len (c.peak_fragments_map[j]) == 0):
                    count += 1
            if (len(c.peaks) > 0):
                self.assertFalse(count == len(c.peaks))

            # test if the peak_fragments_map is correct and complete
            if sample['type'] == 'synthetic':
                self.assertTrue(count == 0)
                mol = Chem.MolFromSmiles(compound['Smiles'])
                fragments = FragmentEngine(Chem.MolToMolBlock(mol), 2, 2, 1, 0, 0)
                fragments.generate_fragments()
                for i, peak in enumerate(c.peaks):
                    annotations = fragments.find_fragments(peak[0], 0.1, self.args['ppm'], self.args['mz_tolerance'])
                    annotations = set([annotation[0] for annotation in annotations])
                    for annotation in annotations:
                        self.assertTrue(annotation in c.peak_fragments_map[i])
                    self.assertTrue(len(c.peak_fragments_map[i]) == len(annotations))
    
    # def test_apply_sirius(self):
    #     # print count of files in helperDirectory
    #     file_not_found = 0
    #     test_count = 100
    #     test_done = 0
    #     for compound in self.data_dict_filtered:
    #         if test_done == test_count:
    #             break
    #         try:
    #             c = Compound(self.data_dict_filtered[compound], args=self.args)
    #         except ValueError as e:
    #             # check if the error is due to an unsupported adduct
    #             if (str(e).startswith("Adduct not supported:")):
    #                 continue
    #             else:
    #                 raise e
    #         if (c.Adduct != "M+H"):
    #             continue
            
    #         before = copy.deepcopy(c.peak_fragments_map)
    #         try:
    #             with open(os.path.join(self.helperDirectory, compound + "_fragmentationtree.json")) as f:
    #                 molSirius = json.load(f)
    #         except FileNotFoundError:
    #             file_not_found += 1
    #             continue
    #         c.apply_sirius(molSirius)
    #         for i, peak in enumerate(c.peaks):
    #             self.assertTrue(len(c.peak_fragments_map[i]) <= len(before[i]))
    #             if (len(c.peak_fragments_map[i]) == 0):
    #                 self.assertTrue(len(before[i]) == 0)
            
    #         test_done += 1
    #     if (file_not_found > test_count/10):
    #         print("WARNING: " + str(file_not_found) + " files were not found in " + self.helperDirectory)

if __name__ == '__main__':
    unittest.main()