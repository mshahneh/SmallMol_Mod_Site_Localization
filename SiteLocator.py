import utils as utils
import visualizer
import numpy as np
from rdkit import Chem
import fragmentation_py as fragmentation_py
from alignment import _cosine_fast


class SiteLocator():
    def __init__(self, molData, modifData, mol, args = {}):

        # set default arguments
        self.args = {
            'adduct': 'M+H',
            'filter_peaks_method': 'top_k',
            'filter_peaks_variable': 50,
            'mz_tolerance': 0.05,
            'ppm': 1.01,
            'min_score_ratio': 0.5,
            'distance_decay': 0.1,
        }
        self.args.update(args)

        if type(molData) != type(modifData):
            raise Exception('molData and modifData must be of the same type')
        if type(molData) == str:
            molUsi = molData
            modifUsi = modifData
            alignment = utils.getMatchedPeaks(molUsi, modifUsi)
            molData = alignment['spectrum1']
            modifData = alignment['spectrum2']
            

        self.molPrecursorMz = molData['precursor_mz']
        self.molPrecursorCharge = molData['precursor_charge']
        self.modifPrecursorMz = modifData['precursor_mz']
        self.modifPrecursorCharge = modifData['precursor_charge']
        self.molPeaks = utils.filter_peaks(molData['peaks'], self.args['filter_peaks_method'], self.args['filter_peaks_variable'])
        self.modifPeaks = utils.filter_peaks(modifData['peaks'], self.args['filter_peaks_method'], self.args['filter_peaks_variable'])
        
        cosine, matchedPeaks = _cosine_fast(utils.convert_to_SpectrumTuple(self.molPeaks, self.molPrecursorMz, self.molPrecursorCharge), 
                                            utils.convert_to_SpectrumTuple(self.modifPeaks, self.modifPrecursorMz, self.modifPrecursorCharge),
                                            self.args['mz_tolerance'], True)
        try:
            self.cosine = alignment['cosine']
        except:
            self.cosine = cosine
        self.matchedPeaks = matchedPeaks

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
    
    def calculate_score(self, peak_presence_only = False, consider_intensity = False):
        self.appearance_shifted = {i: [] for i in range(0, self.molMol.GetNumAtoms())}
        self.appearance_unshifted = {i: [] for i in range(0, self.molMol.GetNumAtoms())}

        shifted = self.shifted
        unshifted = self.unshifted

        scores_unshifted = np.zeros(self.molMol.GetNumAtoms())
        scores_shifted = np.zeros(self.molMol.GetNumAtoms())
        max_intensity = max([_[1] for _ in self.molPeaks])

        for peak in unshifted:
            atom_presence = set()
            possiblities = self.fragments.find_fragments(self.molPeaks[peak[0]][0], 0.1, self.args['ppm'], self.args['mz_tolerance'])
            for possibility in possiblities:
                smiles = self.fragments.get_fragment_info(possibility[0], 0)[3]
                substructure = Chem.MolFromSmiles(smiles, sanitize=False)
                hitAtoms, hitBonds = utils.getHitAtomsAndBonds(self.molMol, substructure)
                for atomSet in hitAtoms:
                    for atom in atomSet:
                        atom_presence.add(atom)
                        self.appearance_unshifted[atom].append(possibility[0])
                        
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
            possiblities = self.fragments.find_fragments(self.molPeaks[peak[0]][0], 0.1, self.args['ppm'], self.args['mz_tolerance'])
            for possibility in possiblities:
                smiles = self.fragments.get_fragment_info(possibility[0], 0)[3]
                substructure = Chem.MolFromSmiles(smiles, sanitize=False)
                hitAtoms, hitBonds = utils.getHitAtomsAndBonds(self.molMol, substructure)
                for atomSet in hitAtoms:
                    for atom in atomSet:
                        atom_presence.add(atom)
                        self.appearance_shifted[atom].append(possibility[0])

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
    
    def accuracy_score(self, modificationSiteIdx, peak_presence_only = False, combine = False, return_all = False, consider_intensity = False):
        
        scores_unshifted, scores_shifted = self.calculate_score(peak_presence_only, consider_intensity)
        scores = self.distance_score( scores_unshifted, scores_shifted, combine)

        maxScore = max(scores)
        if maxScore == 0:
            if return_all:
                return {'score': 0, 'count': 0, 'isMax': 0, 'closestMaxAtomDistance': 0}
            else:
                return 0
        
        for i in range(self.molMol.GetNumAtoms()):
            if scores[i] < self.args['min_score_ratio'] * maxScore:
                scores[i] = 0

        scores /= np.sum(scores)
        maxScore = max(scores)
        
        closestMaxAtomIndx = 0
        localDistances = 0
        count = 0
        graphDiameter = np.amax(self.distances)


        entropy = 0
        
        weights_distance = np.exp(-self.distances[modificationSiteIdx] * self.args['distance_decay'])
        scores_wegihted = np.dot(scores, weights_distance)
        H_weighted = -np.sum(np.dot(scores_wegihted, np.log(scores)))


        

        for i in range(0, self.molMol.GetNumAtoms()):
            if scores[i] > 0:
                count += 1
                localDistances += (self.distances[modificationSiteIdx][i]/graphDiameter) * scores[i]/maxScore
                if scores[i] == maxScore and self.distances[modificationSiteIdx][i] < self.distances[modificationSiteIdx][closestMaxAtomIndx]:
                    closestMaxAtomIndx = i
        
        if return_all:
            res = {'score': 1 - (localDistances/count)}
            res['count'] = count
            res['isMax'] = 1 if scores[modificationSiteIdx] == maxScore else 0
            res['closestMaxAtomDistance'] = self.distances[modificationSiteIdx][closestMaxAtomIndx]
        else:
            res = 1 - (localDistances/count)
        
        return res
    
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
        structures = []
        possiblities = self.fragments.find_fragments(peak_weight, 0.1, 1, mz_precision_abs)
        for possibility in possiblities:
            smiles = self.fragments.get_fragment_info(possibility[0], 0)[3]
            substructure = Chem.MolFromSmiles(smiles, sanitize=False)
            if self.molMol.HasSubstructMatch(substructure):
                structures.append(smiles)
            else:
                print("Error handling substructure, no match found: ", smiles, "")
        return structures



