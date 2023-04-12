import utils as utils
import visualizer
import numpy as np
import copy
from rdkit import Chem
import fragmentation_py as fragmentation_py
from alignment import _cosine_fast, SpectrumTuple, handle_alignment


class SiteLocator():
    def __init__(self, molData, modifData, mol, args = {}):

        # set default arguments
        self.args = {
            'adduct': 'M+H',
            'filter_peaks_method': 'top_k',
            'filter_peaks_variable': 50,
            'mz_tolerance': 0.05,
            'ppm': 1.01,
        }
        self.args.update(args)
        print (self.args)

        res = handle_alignment(molData, modifData, self.args)
        molData = res['molData']
        modifData = res['modifData']
        
        
        self.cosine = res['cosine']
        self.matchedPeaks = res['matchedPeaks']
        self.molData = molData
        self.modifData = modifData   

        self.molPrecursorMz = molData['precursor_mz']
        self.molPrecursorCharge = molData['precursor_charge']
        self.modifPrecursorMz = modifData['precursor_mz']
        self.modifPrecursorCharge = modifData['precursor_charge']
        self.molPeaks = molData['peaks']
        self.modifPeaks = modifData['peaks']

        if type(mol) == str:
            self.molMol = Chem.MolFromSmiles(mol)
        else:
            self.molMol = mol
        
        self.distances = Chem.rdmolops.GetDistanceMatrix(self.molMol)
        
        self.fragments = fragmentation_py.FragmentEngine(Chem.MolToMolBlock(self.molMol), 2, 2, 1, 0, 0)
        self.numFrag = self.fragments.generate_fragments()
          
        shifted, unshifted = utils.separateShifted(self.matchedPeaks, self.molPeaks, self.modifPeaks)
        self.shifted = shifted # list of tuples (molPeakIndex, modifPeakIndex) of matched shifted peaks
        self.unshifted = unshifted # list of tuples (molPeakIndex, modifPeakIndex) of matched unshifted peaks
        self.appearance_shifted = None # dictionary of lists of shifted peak fragments for each atom
        self.appearance_unshifted = None # dictionary of lists of unshifted peak fragments for each atom

        self.mol_filtered_annotations = dict()

    def update_filtered_annotations(self, peak, fragments):
        if fragments is None or len(fragments) == 0:
            # remove peak from self.mol_filtered_annotations
            self.mol_filtered_annotations.pop(peak, None)
            return
        
        if peak in self.mol_filtered_annotations:
            # intersection of the two lists
            self.mol_filtered_annotations[peak] = list(set(self.mol_filtered_annotations[peak]) & set(fragments))
        else:
            self.mol_filtered_annotations[peak] = fragments
    
    def helper_molecule(self, helperMolData, helperMolStruct):

        res = handle_alignment(self.molData, helperMolData, self.args)
        helperMolData = res['modifData']
        if type(helperMolStruct) == str:
            helperMolStruct = Chem.MolFromSmiles(helperMolStruct)
        
        if not self.molMol.HasSubstructMatch(helperMolStruct) or not helperMolStruct.HasSubstructMatch(self.molMol):
            print("helper molecule is not a substructure of the main molecule")
            return

        matechedPeaks = res['matchedPeaks']
        shifted, unshifted = utils.separateShifted(matechedPeaks, self.molPeaks, helperMolData['peaks'])
        helperFrag = fragmentation_py.FragmentEngine(Chem.MolToMolBlock(helperMolStruct), 2, 2, 1, 0, 0)
        helperFrag.generate_fragments()

        for peak in unshifted:
            ## update self.mol_filtered_annotations
            # only annotations that are in both molecules
            molAnnotations = self.fragments.find_fragments(self.molPeaks[peak[0]][0], 0.1, self.args['ppm'], self.args['mz_tolerance'])
            helperAnnotations = helperFrag.find_fragments(helperMolData['peaks'][peak[1]][0], 0.1, self.args['ppm'], self.args['mz_tolerance'])
            posibilities = set()
            flag = False
            for molAnnotation in molAnnotations:
                for helperAnnotation in helperAnnotations:
                    molSubSmiles = self.fragments.get_fragment_info(molAnnotation[0], 0)[3]
                    helperSubSmiles = helperFrag.get_fragment_info(helperAnnotation[0], 0)[3]
                    if Chem.CanonSmiles(molSubSmiles) == Chem.CanonSmiles(helperSubSmiles):
                        posibilities.add(molAnnotation[0])
                        flag = True
                        break
                if not flag:
                    print("there wasn't any match for this peak", peak, molAnnotation)
            
            # print ("unshi", peak, posibilities, molAnnotations, helperAnnotations)
            self.update_filtered_annotations(peak[0], posibilities)

        for peak in shifted:
            ## update self.mol_filtered_annotations
            # check annotations for each peak to see if they contain the modification site
            molAnnotations = self.fragments.find_fragments(self.molPeaks[peak[0]][0], 0.1, self.args['ppm'], self.args['mz_tolerance'])
            helperAnnotations = helperFrag.find_fragments(helperMolData['peaks'][peak[1]][0], 0.1, self.args['ppm'], self.args['mz_tolerance'])
            posibilities = set()
            for molAnnotation in molAnnotations:
                for helperAnnotation in helperAnnotations:
                    molSubSmiles = self.fragments.get_fragment_info(molAnnotation[0], 0)[3]
                    helperSubSmiles = helperFrag.get_fragment_info(helperAnnotation[0], 0)[3]
                    molSubStruct = Chem.MolFromSmiles(molSubSmiles)
                    helperSubStruct = Chem.MolFromSmiles(helperSubSmiles)
                    if self.molPrecursorMz < helperMolData['precursor_mz']:
                        if helperSubStruct.HasSubstructMatch(molSubStruct):# and not self.molMol.HasSubstructMatch(helperSubStruct):
                            posibilities.add(molAnnotation[0])
                            break
                    else:
                        if molSubStruct.HasSubstructMatch(helperSubStruct):# and not helperMolStruct.HasSubstructMatch(molSubStruct):
                            posibilities.add(molAnnotation[0])
                            break
            
            # print ("shifted", peak, posibilities, res['modifData']['peaks'][peak[0]][0], helperMolData['peaks'][peak[1]][0])
            self.update_filtered_annotations(peak[0], posibilities)

    
    def calculate_score(self, peak_presence_only = False, consider_intensity = False, modificationSite = None):
        self.appearance_shifted = {i: [] for i in range(0, self.molMol.GetNumAtoms())}
        self.appearance_unshifted = {i: [] for i in range(0, self.molMol.GetNumAtoms())}

        shifted = self.shifted
        unshifted = self.unshifted

        scores_unshifted = np.zeros(self.molMol.GetNumAtoms())
        scores_shifted = np.zeros(self.molMol.GetNumAtoms())
        max_intensity = max([_[1] for _ in self.molPeaks])

        for peak in unshifted:
            atom_presence = set()
            if peak[0] in self.mol_filtered_annotations:
                possiblities = list(self.mol_filtered_annotations[peak[0]])
            else:
                possiblities = self.fragments.find_fragments(self.molPeaks[peak[0]][0], 0.1, self.args['ppm'], self.args['mz_tolerance'])
                possiblities = set([x[0] for x in possiblities])
            
            for possibility in possiblities:
                hitAtoms = self.fragments.get_fragment_info(possibility, 0)[1]
                for atom in hitAtoms:
                    atom_presence.add(atom)
                    self.appearance_unshifted[atom].append(possibility)
                    
                    if not peak_presence_only:
                        if consider_intensity:
                            scores_unshifted[atom] += self.molPeaks[peak[0]][1]/max_intensity
                        else:
                            scores_unshifted[atom] += 1
            
            if peak_presence_only:
                for atom in atom_presence:
                    if consider_intensity:
                        scores_unshifted[atom] += self.molPeaks[peak[0]][1]/max_intensity
                    else:
                        scores_unshifted[atom] += 1
            
        for peak in shifted:
            atom_presence = set()
            if peak[0] in self.mol_filtered_annotations:
                possiblities = list(self.mol_filtered_annotations[peak[0]])
            else:
                possiblities = self.fragments.find_fragments(self.molPeaks[peak[0]][0], 0.1, self.args['ppm'], self.args['mz_tolerance'])
                possiblities = set([x[0] for x in possiblities])

            for possibility in possiblities:
                smiles = self.fragments.get_fragment_info(possibility, 0)[3]
                substructure = Chem.MolFromSmiles(smiles, sanitize=False)
                hitAtoms = self.fragments.get_fragment_info(possibility, 0)[1]
                if modificationSite is None or modificationSite in hitAtoms:
                    for atom in hitAtoms:
                        atom_presence.add(atom)
                        self.appearance_shifted[atom].append(possibility)

                        if not peak_presence_only:
                            if consider_intensity:
                                scores_shifted[atom] += self.molPeaks[peak[0]][1] / max_intensity
                            else:
                                scores_shifted[atom] += 1
            
            if peak_presence_only:
                for atom in atom_presence:
                    if consider_intensity:
                        scores_shifted[atom] += self.molPeaks[peak[0]][1] / max_intensity
                    else:
                        scores_shifted[atom] += 1

        return scores_unshifted, scores_shifted
    
    def distance_score(self, scores_unshifted, scores_shifted, combine = True):
        scores = np.zeros(self.molMol.GetNumAtoms())
        for i in range(0, self.molMol.GetNumAtoms()):
            if combine:
                scores[i] = scores_shifted[i] - scores_unshifted[i]
            else:
                scores[i] = scores_shifted[i]
        
        ## normalize
        scores = (scores - np.amin(scores)) / (np.amax(scores) - np.amin(scores))

        return scores
    

    def tempScore(self, modificationSiteIdx, preds, return_all = False):
        
        maxScore = max(preds)
        if maxScore == 0:
            if return_all:
                return {'score': 0, 'count': 0, 'isMax': 0, 'closestMaxAtomDistance': 0}
            else:
                return 0
    

        for i in range(self.molMol.GetNumAtoms()):
            if preds[i] < 0/5 * maxScore:
                preds[i] = 0
        preds /= np.sum(preds)
        maxScore = max(preds)
        graphDiameter = np.amax(self.distances)
        count = 0
        localDistances = 0
        closestMaxAtomIndx = 0
        # print("DUAAAM", graphDiameter, self.molMol.GetNumAtoms())
        for i in range(self.molMol.GetNumAtoms()):
            if preds[i] > 0:
                # print("in if")
                count += 1

                # print("ASD", self.distances[modificationSiteIdx][i])
                localDistances += (self.distances[modificationSiteIdx][i]/graphDiameter) * preds[i]/maxScore
                if preds[i] == maxScore and self.distances[modificationSiteIdx][i] < self.distances[modificationSiteIdx][closestMaxAtomIndx]:
                    closestMaxAtomIndx = i
        
        if return_all:
            res = {'score': 1 - (localDistances/count)}
            res['count'] = count
            res['isMax'] = 1 if preds[modificationSiteIdx] == maxScore else 0
            res['closestMaxAtomDistance'] = self.distances[modificationSiteIdx][closestMaxAtomIndx]
        else:
            res = 1 - (localDistances/count)
        
        return res

    def accuracy_score(self, modificationSiteIdx, peak_presence_only = False, combine = False, return_all = False, consider_intensity = False):
        scores_unshifted, scores_shifted = self.calculate_score(peak_presence_only, consider_intensity)
        scores = self.distance_score(scores_unshifted, scores_shifted, combine)
        return self.tempScore(modificationSiteIdx, scores, return_all)
    
    def get_max_possible_score(self, modificationSiteIdx, peak_presence_only = False, combine = False, return_all = False, consider_intensity = False):
        scores_unshifted, scores_shifted = self.calculate_score(peak_presence_only, consider_intensity, modificationSite=modificationSiteIdx)
        scores = self.distance_score(scores_unshifted, scores_shifted, combine)
        return self.tempScore(modificationSiteIdx, scores, return_all)
    
    def get_peaks_per_atom(self, atomIdx):
        return self.appearance_shifted[atomIdx], self.appearance_unshifted[atomIdx]

    def get_structures_per_atom(self, atomIdx):
        structures_shifted = []
        structures_unshifted = []
        for frag in self.appearance_shifted[atomIdx]:
            structures_shifted.append(self.get_fragment_info(frag, 0)[3])
        for frag in self.appearance_unshifted[atomIdx]:
            structures_unshifted.append(self.get_fragment_info(frag, 0)[3])
        return structures_shifted, structures_unshifted

    def get_structures_per_peak(self, peak_weight, mz_precision_abs = 0.05):
        # start debugging
        print ("debugging", peak_weight)
        possiblities = self.fragments.find_fragments(peak_weight, 0.1, self.args['ppm'], self.args['mz_tolerance'])
        for possibility in possiblities:
            print("debugging", possibility, self.fragments.get_fragment_info(possibility[0], 0))

        # end debugging
        structures = []
        result_posibility_indicies = []
        ind = -1
        for i in range(len(self.molPeaks)):
            if abs(round(self.molPeaks[i][0], 4) - peak_weight) < 0.00001:
                ind = i
        print("ind", ind)
        print(self.mol_filtered_annotations)
        if ind in self.mol_filtered_annotations:
            print("filtered", ind, self.mol_filtered_annotations[ind])
            possiblities = list(self.mol_filtered_annotations[ind])
        else:
            possiblities = self.fragments.find_fragments(peak_weight, 0.1, self.args['ppm'], self.args['mz_tolerance'])
            possiblities = set([x[0] for x in possiblities])
        print("get_structures_per_peak", possiblities)
        for possibility in possiblities:
            print(self.fragments.get_fragment_info(possibility, 0))
            smiles = self.fragments.get_fragment_info(possibility, 0)[3]
            posibility_indices = self.fragments.get_fragment_info(possibility, 0)[1]
            substructure = Chem.MolFromSmiles(smiles, sanitize=False)
            if self.molMol.HasSubstructMatch(substructure):
                structures.append(smiles)
                result_posibility_indicies.append(posibility_indices)
            else:
                print("Error handling substructure, no match found: ", smiles, "")
        return structures, result_posibility_indicies



