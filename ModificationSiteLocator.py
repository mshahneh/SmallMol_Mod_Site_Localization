from .alignment_n import align
import json
from . import utils_n as utils
import numpy as np
from . import calculate_scores_n as Calc_Scores
from . import Compound_n as Compound
from rdkit import Chem

class ModificationSiteLocator():
    def __init__(self, main_compound, modified_compound, args = {}):
        
        self.args = {"mz_tolerance": 0.1,
                      "ppm": 40,
                    }
        self.args.update(args)

        for arg in list(args.keys()):
            if type(args[arg]) == str:
                try:
                    self.args[arg] = float(args[arg])
                except:
                    self.args[arg] = args[arg]
        
        if type(main_compound) == str:
            main_compound = Compound.Compound(main_compound, args=self.args)
        
        if type(modified_compound) == str:
            modified_compound = Compound.Compound(modified_compound, args=self.args)

        self.main_compound = main_compound
        self.modified_compound = modified_compound
        self.main_compound.remove_large_peaks()
        self.modified_compound.remove_large_peaks()
        self.cosine, self.matched_peaks = align(self.main_compound, self.modified_compound, self.args["mz_tolerance"], self.args["ppm"])
        self.shifted, self.unshifted = utils.separateShifted(self.matched_peaks, 
                                                                 self.main_compound.peaks, self.modified_compound.peaks)

    
    def __repr__(self) -> str:
        obj = {"main_compound": self.main_compound, "modified_compound": self.modified_compound}
        return json.dumps(obj)
    
    def find_existance(self, peakids):
        """
        For the given peaks, finds the fragments that each atom contributes to
        input:
            peakids: list of peak ids
        output:
            existance: list of dicts for each atom, each dict contains the peak ids as keys and the fragments as values
        """
        existance = [dict() for i in range(len(self.main_compound.structure.GetAtoms()))]
        for peak in peakids:
            for fragment in self.main_compound.peak_fragments_map[peak]:
                hitAtoms = self.main_compound.fragments.get_fragment_info(fragment, 0)[1]
                for atom in hitAtoms:
                    if peak not in existance[atom]:
                        existance[atom][peak] = []
                    existance[atom][peak].append(fragment)
        return existance
    
    def calculate_contributions(self, peakids, CI = False, CPA = True, CFA = True):
        """ 
        input:
            peakids: list of peak ids
            CI: (Consider_Intensity) bool, if True, the intensity of the peaks is considered (default: False)
            CPA: (Consider_Peak_Ambiguity) bool, if True, the peak ambiguity (number of fragments assigned to a peak) is considered (default: True)
            CFA: (Consider_Fragment_Ambiguity) bool, if True, the fragment ambiguity (number of atoms in fragment) is considered (default: True)
        """
        existance_data = self.find_existance(peakids)
        contributions = [0 for i in range(len(self.main_compound.structure.GetAtoms()))]
        for atom in range(len(existance_data)):
            for peak in existance_data[atom]:
                intensity_factor = 1
                atom_peak_ambiguity_factor = 1
                fragment_ambiguity_factor = 1

                if CI:
                    intensity_factor = self.main_compound.peaks[peak][1]
                if CPA:
                    atom_peak_ambiguity_factor = len(existance_data[atom][peak])/len(self.main_compound.peak_fragments_map[peak])
                
                for frag in existance_data[atom][peak]:
                    if CFA:
                        fragment_ambiguity_factor = len(self.main_compound.fragments.get_fragment_info(frag, 0)[1])/len(self.main_compound.structure.GetAtoms())
                    
                    contributions[atom] += intensity_factor * atom_peak_ambiguity_factor * fragment_ambiguity_factor
        
        return contributions

    
    def generate_probabilities(self, shifted_only = True, CI = False, CPA = True, CFA = True, method = "old"):
        """"Generate the probabilities for each atom to be the modification site.
        input:
            shifted_only: bool, if True, only the shifted peaks are considered
            CI: (Consider_Intensity) bool, if True, the intensity of the peaks is considered (default: False)
            CPA: (Consider_Peak_Ambiguity) bool, if True, the peak ambiguity (number of fragments assigned to a peak) is considered when calculating the contribution (default: True)
            CFA: (Consider_Fragment_Ambiguity) bool, if True, the fragment ambiguity (number of atoms in fragment) is considered (default: True)
        """
        
        if method == "random_choice":
            n = len(self.main_compound.structure.GetAtoms())
            probabilities = np.zeros(n)
            random_choice = np.random.choice(n)
            probabilities[random_choice] = 1
            return probabilities
        elif method == "random_distribution":
            probabilities = np.random.rand(len(self.main_compound.structure.GetAtoms()))
            probabilities = probabilities / np.sum(probabilities)
            return probabilities
        elif method == "all_equal":
            probabilities = np.ones(len(self.main_compound.structure.GetAtoms()))
            probabilities = probabilities / np.sum(probabilities)
            return probabilities
        elif method == "random_skewed":
            target =  np.random.choice(len(self.main_compound.structure.GetAtoms()))
            probabilities = np.random.rand(len(self.main_compound.structure.GetAtoms()))
            probabilities = probabilities * self.main_compound.distances[target]
            probabilities = probabilities - np.min(probabilities)
            probabilities = probabilities / np.sum(probabilities)
            return probabilities
        elif method.startswith("multiple"):
            # set seed to be 0
            np.random.seed(0)
            res = []
            method_name = method[9:]
            for i in range(50):
                res.append(self.generate_probabilities(method=method_name))
            return res
        elif method == "multiply":
            probabilities = np.zeros(len(self.main_compound.structure.GetAtoms()))
            for atom in range(len(self.main_compound.structure.GetAtoms())):
                positive_contributions = 1
                for peak in self.shifted:
                    count = 0
                    for fragment in self.main_compound.peak_fragments_map[peak[0]]:
                        if atom in self.main_compound.fragments.get_fragment_info(fragment, 0)[1]:
                            count += 1
                    eps = 0
                    if CPA:
                        eps = 1
                    
                    if (len(self.main_compound.peak_fragments_map[peak[0]]) + eps) == 0:
                        positive_contributions *= 0
                    else:
                        positive_contributions *= (count + eps)/(len(self.main_compound.peak_fragments_map[peak[0]]) + eps)
                

                if shifted_only:
                    negative_contributions = 0
                else:
                    # TODO: implement this
                    negative_contributions = 1
                
                probabilities[atom] = positive_contributions - negative_contributions
            
        else:
            s_peakids = [_[0] for _ in self.shifted]
            positive_contributions = self.calculate_contributions(s_peakids, CI=CI, CPA=CPA, CFA=CFA)
            if not shifted_only:
                u_peakids = [_[0] for _ in self.unshifted]
                negative_contributions = self.calculate_contributions(u_peakids, CI=CI, CPA=CPA, CFA=CFA)
            else:
                negative_contributions = [0 for i in range(len(self.main_compound.structure.GetAtoms()))]
            
            probabilities = np.zeros(len(self.main_compound.structure.GetAtoms()))
            for i in range(len(positive_contributions)):
                probabilities[i] = positive_contributions[i] - negative_contributions[i]
            probabilities = Calc_Scores.power_prob(probabilities) 

        return probabilities
    
    def calculate_score(self, true_modification_site, method, probabilities = None, extensive_response = False, filter_ratio = 0.5):
        """Calculate the score for a probability.
        input:
            true_modification_site: int, the true modification site
            method: str, the method to calculate the score
            probabilities: list, the probabilities for each atom to be the modification site
            extensive_response: bool, if True, the function returns a dictionary with more information
            filter_ratio: float, the ratio of the highest probability that is considered 
                        (if 0, all probabilities are considered, if 1, only the highest probability is considered)
        """
        if probabilities is None:
            probabilities = self.generate_probabilities()
        
        G = self.main_compound.distances
        return Calc_Scores.calculate(G, probabilities, true_modification_site, method)
    
    def get_structures_by_peak_index(self, peakindex):
        """Get all the annotations for a peak."""
        structures = []
        structure_indicies = []
        frags = []
        for fragment in self.main_compound.peak_fragments_map[peakindex]:
            fragInfo = self.main_compound.fragments.get_fragment_info(fragment, 0)
            smiles = fragInfo[3]
            hitAtoms = fragInfo[1]
            structures.append(smiles)
            structure_indicies.append(hitAtoms)
            frags.append(fragment)
        
        return structures, structure_indicies, frags

    
    def get_structures_by_peak_weight(self, peak_weight, mz_precision_abs = None, mz_ppm_tolerance = None):
        """Get all the annotations for a peak weight."""
        if mz_precision_abs is None:
            mz_precision_abs = self.args['mz_tolerance']
        if mz_ppm_tolerance is None:
            mz_ppm_tolerance = self.args['ppm']
        structures = []
        structure_indicies = []
        fragments = []
        ind = self.main_compound.get_peak_index(peak_weight, {"mz_tolerance": mz_precision_abs, "ppm": mz_ppm_tolerance})
        for i in ind:
            structures_i, structure_indicies_i, fragments_i = self.get_structures_by_peak_index(i)
            structures += structures_i
            structure_indicies += structure_indicies_i
            fragments += fragments_i
        return structures, structure_indicies, fragments
    
    def get_main_shifted_index(self):
        """Get the indexes of the shifted peaks in the main compound."""
        return [_[0] for _ in self.shifted]
    
    def get_main_unshifted_index(self):
        """Get the indexes of the unshifted peaks in the main compound."""
        return [_[0] for _ in self.unshifted]
    
    def get_modified_shifted_index(self):
        """Get the index of the shifted peaks in the modified compound."""
        return [_[1] for _ in self.shifted]
    
    def get_modified_unshifted_index(self):
        """Get the index of the unshifted peaks in the modified compound."""
        return [_[1] for _ in self.unshifted]
    
    def apply_helpers(self, helpers, unshifted_mode = None):
        """Apply the helpers to the main compound."""
        for helper in helpers:
            helper_compound = Compound.Compound(helper, args=self.args)
            cosine, matched_peaks = align(self.main_compound, helper_compound, self.args["mz_tolerance"], self.args["ppm"])
            shifted, unshifted = utils.separateShifted(matched_peaks, self.main_compound.peaks, helper_compound.peaks)
            self.main_compound.apply_helper(helper_compound, shifted, unshifted, unshifted_mode)
