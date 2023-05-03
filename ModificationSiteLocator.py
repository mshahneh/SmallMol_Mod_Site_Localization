from alignment_n import align
import json
import utils as utils_old
import numpy as np
import calculate_scores_n as Calc_Scores

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
    
    def calculate_contributions(self, peakids, PPO = False, CO = False, true_modification_site = None):
        """ 
        input:
            peakids: list of peak ids
            PPO: (Peak_Presense_only) bool, if True, only the number of contributed peaks is considered, not the number of fragments
            CO: (Consider_Intensity) bool, if True, the intensity of the peaks is considered (default: False)
            true_modification_site: int, if not None, only fragments that contain this atom are considered
        """
        contributions_data = self.find_contributions(peakids, true_modification_site)
        contributions = [0 for i in range(len(self.main_compound.structure.GetAtoms()))]
        intensity_factor = 0
        if CO:
            intensity_factor = 1
        for atom in range(len(contributions_data)):
            for peak in contributions_data[atom]:
                if PPO:
                    contributions[atom] += 1 * intensity_factor * self.main_compound.peaks[peak][1]
                else:
                    contributions[atom] += len(contributions_data[atom][peak]) * intensity_factor * self.main_compound.peaks[peak][1]
        return contributions

    
    def generate_probabilities(self, shifted_only = False, PPO = False, CO = False, true_modification_site = None):
        """"Generate the probabilities for each atom to be the modification site.
        input:
            shifted_only: bool, if True, only the shifted peaks are considered
            PPO: (Peak_Presense_only) bool, if True, only the number of contributed peaks is considered, not the number of fragments
            CO: (Consider_Intensity) bool, if True, the intensity of the peaks is considered (default: False)
        """
        s_peakids = [_[0] for _ in self.shifted]
        positive_contributions = self.calculate_contributions(s_peakids, PPO, CO, true_modification_site)
        if not shifted_only:
            u_peakids = [_[0] for _ in self.shifted]
            negative_contributions = self.calculate_contributions(u_peakids, PPO, CO, None)
        else:
            negative_contributions = [0 for i in range(len(self.main_compound.structure.GetAtoms()))]
        
        probabilities = np.zeros(len(self.main_compound.structure.GetAtoms()))
        for i in range(len(positive_contributions)):
            probabilities[i] = positive_contributions[i] - negative_contributions[i]
        
        # Normalize probabilities
        probabilities = probabilities - np.min(probabilities)
        probabilities = probabilities / np.sum(probabilities)
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
        if maxScore == 0:
            if extensive_response:
                return {'score': 0, 'count': 0, 'isMax': 0, 'closestMaxAtomDistance': 0}
            else:
                return 0
        
        for i in range(self.main_compound.structure.GetNumAtoms()):
            if probabilities[i] <= filter_ratio * maxScore:
                probabilities[i] = 0
        probabilities /= np.sum(probabilities)
        
        # call the score function from Calc_Scores based on the method
        if method == "is_max":
            return Calc_Scores.is_max(G, probabilities, true_modification_site)
        elif method == "dist_from_max":
            return Calc_Scores.dist_from_max(G, probabilities, true_modification_site)
        elif method == "average_dist_from_max":
            return Calc_Scores.average_dist_from_max(G, probabilities, true_modification_site)
        elif method == "average_dist":
            return Calc_Scores.average_dist(G, probabilities, true_modification_site)
        else:
            raise Exception("Method not found")
        
    
    def get_structures_per_peak(self, peak_weight, mz_precision_abs = 0.05):
        pass
        # structures = []
        # result_posibility_indicies = []
        # ind = -1
        # for i in range(len(self.molPeaks)):
        #     if abs(round(self.molPeaks[i][0], 4) - peak_weight) < 0.00001:
        #         ind = i
        # if ind in self.mol_filtered_annotations:
        #     possiblities = list(self.mol_filtered_annotations[ind])
        # else:
        #     possiblities = self.fragments.find_fragments(peak_weight, 0.1, self.args['ppm'], self.args['mz_tolerance'])
        #     possiblities = set([x[0] for x in possiblities])
        # for possibility in possiblities:
        #     smiles = self.fragments.get_fragment_info(possibility, 0)[3]
        #     posibility_indices = self.fragments.get_fragment_info(possibility, 0)[1]
        #     substructure = Chem.MolFromSmiles(smiles, sanitize=False)
        #     if self.molMol.HasSubstructMatch(substructure):
        #         structures.append(smiles)
        #         result_posibility_indicies.append(posibility_indices)
        #     else:
        #         if self.verbose > 0:
        #             print("Error handling substructure, no match found: ", smiles, "")
        # return structures, result_posibility_indicies