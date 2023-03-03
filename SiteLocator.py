import utils as utils
import visualizer
import numpy as np
from rdkit import Chem
import fragmentation_py as fragmentation_py


class SiteLocator():
    def __init__(self, molUsi, modifUsi, mol):
        self.molUsi = molUsi #smaler molecule
        self.modifUsi = modifUsi #bigger (with addition) modified molecule
        if type(mol) == str:
            self.molMol = Chem.MolFromSmiles(mol, sanitize=False)
        else:
            self.molMol = mol
        alignment = utils.getMatchedPeaks(molUsi, modifUsi)
        self.molMeta = alignment['spectrum1']
        self.molPeaks = alignment['spectrum1']['peaks']
        self.modifMeta = alignment['spectrum2']
        self.modifPeaks = alignment['spectrum2']['peaks']
        self.cosine = alignment['cosine']
        self.matchedPeaks = alignment['peak_matches']
        
        self.distances = Chem.rdmolops.GetDistanceMatrix(self.molMol)
        
        self.appearance_shifted = {i: [] for i in range(0, self.molMol.GetNumAtoms())}
        self.appearance_unshifted = {i: [] for i in range(0, self.molMol.GetNumAtoms())}
        self.fragments = fragmentation_py.FragmentEngine(Chem.MolToMolBlock(self.molMol), 3, 2, 1, 0, 0)
          
        shifted, unshifted = utils.separateShifted(self.matchedPeaks, self.molPeaks, self.modifPeaks)
        self.shifted = shifted
        self.unshifted = unshifted
    
    def calculate_score(self):
        self.appearance_shifted = {i: [] for i in range(0, self.molMol.GetNumAtoms())}
        self.appearance_unshifted = {i: [] for i in range(0, self.molMol.GetNumAtoms())}

        shifted, unshifted = utils.separateShifted(self.matchedPeaks, self.molPeaks, self.modifPeaks)
        
        
        numFrag = self.fragments.generate_fragments()
        # print("Number of fragments: ", numFrag)

        for peak in unshifted:
            possiblities = self.fragments.find_fragments(self.molPeaks[peak[0]][0], 0.1, 0.98, 0.5)
            for possibility in possiblities:
                smiles = self.fragments.get_fragment_info(possibility[0], 0)[3]
                substructure = Chem.MolFromSmiles(smiles, sanitize=False)
                hitAtoms, hitBonds = utils.getHitAtomsAndBonds(self.molMol, substructure)
                for atomSet in hitAtoms:
                    for atom in atomSet:
                        self.appearance_unshifted[atom].append(possibility[0])
            
        for peak in shifted:
            possiblities = self.fragments.find_fragments(self.molPeaks[peak[0]][0], 0.1, 0.98, 0.5)
            for possibility in possiblities:
                smiles = self.fragments.get_fragment_info(possibility[0], 0)[3]
                substructure = Chem.MolFromSmiles(smiles, sanitize=False)
                hitAtoms, hitBonds = utils.getHitAtomsAndBonds(self.molMol, substructure)
                for atomSet in hitAtoms:
                    for atom in atomSet:
                        self.appearance_shifted[atom].append(possibility[0])
        
        scores_unshifted = []
        scores_shifted = []
        for i in range(0, self.molMol.GetNumAtoms()):
            scores_unshifted.append(len(self.appearance_unshifted[i]))
            scores_shifted.append(len(self.appearance_shifted[i]))
        return scores_unshifted, scores_shifted
    
    def distance_score(self, scores_unshifted, scores_shifted):
        scores = {}
        for i in range(0, self.molMol.GetNumAtoms()):
            scores[i] = scores_shifted[i] - scores_unshifted[i]
        
        ## add min of scores to all in scores
        minScore = min(scores.values())
        for i in range(0, self.molMol.GetNumAtoms()):
            scores[i] += abs(minScore)
        return scores
    
    def accuracy_score(self, modificationSiteIdx):
        
        scores_unshifted, scores_shifted = self.calculate_score()
        scores = self.distance_score( scores_unshifted, scores_shifted)

        Max = max(scores.values())
        if Max == 0:
            return 0
        if scores[modificationSiteIdx] == Max:
            return 1
        else:
            distances = 0
            count = 0
            for i in range(0, self.molMol.GetNumAtoms()):
                if scores[i] == Max:
                    count += 1
                    distances += (self.distances[modificationSiteIdx][i]/np.amax(self.distances))
            return 1 - (distances / count)
    
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

    def get_structures_per_peak(self, peak_weight):
        structures = []
        possiblities = self.fragments.find_fragments(peak_weight, 0.1, 0.98, 0.5)
        for possibility in possiblities:
            structures.append(self.fragments.get_fragment_info(possibility[0], 0)[3])
        return structures



