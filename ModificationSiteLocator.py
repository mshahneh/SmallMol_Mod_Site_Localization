from alignment_n import align
import json
import utils as utils_old
import numpy as np
import calculate_scores_n as Calc_Scores
from rdkit import Chem

class ModificationSiteLocator():
    def __init__(self, main_compound, modified_compound, args = {}):
        
        self.args = {"mz_tolerance": 0.1}
        self.args.update(args)

        self.main_compound = main_compound
        self.modified_compound = modified_compound
        self.cosine, self.matched_peaks = align(self.main_compound, self.modified_compound, self.args["mz_tolerance"])
        self.shifted, self.unshifted = utils_old.separateShifted(self.matched_peaks, 
                                                                 self.main_compound.peaks, self.modified_compound.peaks)

    
    def __repr__(self) -> str:
        obj = {"main_compound": self.main_compound, "modified_compound": self.modified_compound}
        return json.dumps(obj)
    
    def find_contributions(self, peakids, true_modification_site = None):
        """
        For the given peaks, finds the fragments that each atom contributes to
        input:
            peakids: list of peak ids
            true_modification_site: int, if not None, only fragments that contain this atom are considered
        output:
            contributions: list of dicts for each atom, each dict contains the peak ids as keys and the fragments as values
        """
        contributions = [dict() for i in range(len(self.main_compound.structure.GetAtoms()))]
        for peak in peakids:
            for fragment in self.main_compound.peak_fragments_map[peak]:
                hitAtoms = self.main_compound.fragments.get_fragment_info(fragment, 0)[1]
                if true_modification_site != None and true_modification_site not in hitAtoms:
                    continue
                for atom in hitAtoms:
                    if peak not in contributions[atom]:
                        contributions[atom][peak] = []
                    contributions[atom][peak].append(fragment)
        return contributions
    
    def calculate_contributions(self, peakids, PPO = False, CI = False, true_modification_site = None):
        """ 
        input:
            peakids: list of peak ids
            PPO: (Peak_Presense_only) bool, if True, only the number of contributed peaks is considered, not the number of fragments in each peak
            CI: (Consider_Intensity) bool, if True, the intensity of the peaks is considered (default: False)
            true_modification_site: int, if not None, only fragments that contain this atom are considered
        """
        contributions_data = self.find_contributions(peakids, true_modification_site)
        # print("debugging: contributions_data", contributions_data)
        contributions = [0 for i in range(len(self.main_compound.structure.GetAtoms()))]
        for atom in range(len(contributions_data)):
            for peak in contributions_data[atom]:
                intensity_factor = 1
                if CI:
                    intensity_factor = self.main_compound.peaks[peak][1]
                
                if PPO:
                    contributions[atom] += 1 * intensity_factor
                else:
                    contributions[atom] += len(contributions_data[atom][peak]) * intensity_factor
        return contributions

    
    def generate_probabilities(self, shifted_only = False, PPO = False, CI = False, true_modification_site = None):
        """"Generate the probabilities for each atom to be the modification site.
        input:
            shifted_only: bool, if True, only the shifted peaks are considered
            PPO: (Peak_Presense_only) bool, if True, only the number of contributed peaks is considered, not the number of fragments
            CO: (Consider_Intensity) bool, if True, the intensity of the peaks is considered (default: False)
        """
        s_peakids = [_[0] for _ in self.shifted]
        positive_contributions = self.calculate_contributions(s_peakids, PPO, CI, true_modification_site)
        if not shifted_only:
            u_peakids = [_[0] for _ in self.unshifted]
            negative_contributions = self.calculate_contributions(u_peakids, PPO, CI, None)
        else:
            negative_contributions = [0 for i in range(len(self.main_compound.structure.GetAtoms()))]
        
        probabilities = np.zeros(len(self.main_compound.structure.GetAtoms()))
        for i in range(len(positive_contributions)):
            probabilities[i] = positive_contributions[i] - negative_contributions[i]
        
        # print("debugging: probabilities1", probabilities, positive_contributions, negative_contributions)
        
        # Normalize probabilities
        probabilities = probabilities - np.min(probabilities)
        if np.sum(probabilities) != 0:
            probabilities = probabilities / np.sum(probabilities)
        else:
            probabilities = np.zeros(len(self.main_compound.structure.GetAtoms()))
        
        # print("debugging: probabilities2", probabilities, positive_contributions, negative_contributions)
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
        
        maxScore = max(probabilities)
        # for i in range(self.main_compound.structure.GetNumAtoms()):
        #     if probabilities[i] <= filter_ratio * maxScore:
        #         probabilities[i] = 0
        
        # if np.sum(probabilities) != 0:
        #     probabilities /= np.sum(probabilities)

        maxScore = max(probabilities)  
        if maxScore == 0:
            if extensive_response:
                return {'score': 0, 'count': 0, 'isMax': 0, 'closestMaxAtomDistance': 0}
            else:
                return 0
        
        # call the score function from Calc_Scores based on the method
        if method == "is_max":
            return Calc_Scores.is_max(G, probabilities, true_modification_site)
        elif method == "dist_from_max":
            return Calc_Scores.dist_from_max(G, probabilities, true_modification_site)
        elif method == "average_dist_from_max":
            return Calc_Scores.average_dist_from_max(G, probabilities, true_modification_site)
        elif method == "average_dist":
            return Calc_Scores.average_dist(G, probabilities, true_modification_site)
        elif method == "temp":
            return Calc_Scores.temp_score(G, probabilities, true_modification_site)
        else:
            raise Exception("Method not found")
    
    def get_structures_by_peak_id(self, peakid):
        """Get all the annotations for a peak."""
        structures = []
        structure_indicies = []
        for fragment in self.main_compound.peak_fragments_map[peakid]:
            fragInfo = self.main_compound.fragments.get_fragment_info(fragment, 0)
            smiles = fragInfo[3]
            hitAtoms = fragInfo[1]
            substructure = Chem.MolFromSmiles(smiles, sanitize=False)
            if self.molMol.HasSubstructMatch(substructure):
                structures.append(smiles)
                structure_indicies.append(hitAtoms)
        
        return structures, structure_indicies

    
    def get_structures_by_peak_weight(self, peak_weight, mz_precision_abs = None):
        """Get all the annotations for a peak weight."""
        if mz_precision_abs is None:
            mz_precision_abs = self.args['mz_tolerance']
        structures = []
        structure_indicies = []
        ind = []
        for i in range(len(self.molPeaks)):
            if abs(self.molPeaks[i][0] - peak_weight) < self.args['mz_tolerance']:
                ind.append(i)
        for i in ind:
            structures_i, structure_indicies_i = self.get_structures_by_peak_id(i)
            structures += structures_i
            structure_indicies += structure_indicies_i
        return structures, structure_indicies