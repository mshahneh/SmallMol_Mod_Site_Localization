import utils as utils
import visualizer
import numpy as np
from rdkit import Chem
import fragmentation_py as fragmentation_py


class SiteLocator():
    def __init__(self, molUsi, modifUsi, mol, adduct = 'M+H'):
        self.molUsi = molUsi #smaler molecule
        self.modifUsi = modifUsi #bigger (with addition) modified molecule
        if type(mol) == str:
            self.molMol = Chem.MolFromSmiles(mol)
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
        self.fragments = fragmentation_py.FragmentEngine(Chem.MolToMolBlock(self.molMol), 2, 2, 1, 0, 0)
        self.numFrag = self.fragments.generate_fragments()
          
        shifted, unshifted = utils.separateShifted(self.matchedPeaks, self.molPeaks, self.modifPeaks)
        self.shifted = shifted
        self.unshifted = unshifted
    
    def calculate_score(self, peak_presence_only = False):
        self.appearance_shifted = {i: [] for i in range(0, self.molMol.GetNumAtoms())}
        self.appearance_unshifted = {i: [] for i in range(0, self.molMol.GetNumAtoms())}

        shifted, unshifted = utils.separateShifted(self.matchedPeaks, self.molPeaks, self.modifPeaks)

        atom_peak_presence_shifted = {i: 0 for i in range(0, self.molMol.GetNumAtoms())}
        atom_peak_presence_unshifted = {i: 0 for i in range(0, self.molMol.GetNumAtoms())}

        for peak in unshifted:
            atom_presence = set()
            possiblities = self.fragments.find_fragments(self.molPeaks[peak[0]][0], 0.1, 0.98, 0.01)
            for possibility in possiblities:
                smiles = self.fragments.get_fragment_info(possibility[0], 0)[3]
                substructure = Chem.MolFromSmiles(smiles, sanitize=False)
                hitAtoms, hitBonds = utils.getHitAtomsAndBonds(self.molMol, substructure)
                for atomSet in hitAtoms:
                    for atom in atomSet:
                        atom_presence.add(atom)
                        self.appearance_unshifted[atom].append(possibility[0])
            
            for atom in atom_presence:
                atom_peak_presence_unshifted[atom] += 1
            
        for peak in shifted:
            atom_presence = set()
            possiblities = self.fragments.find_fragments(self.molPeaks[peak[0]][0], 0.1, 0.98, 0.5)
            for possibility in possiblities:
                smiles = self.fragments.get_fragment_info(possibility[0], 0)[3]
                substructure = Chem.MolFromSmiles(smiles, sanitize=False)
                hitAtoms, hitBonds = utils.getHitAtomsAndBonds(self.molMol, substructure)
                for atomSet in hitAtoms:
                    for atom in atomSet:
                        atom_presence.add(atom)
                        self.appearance_shifted[atom].append(possibility[0])
            
            for atom in atom_presence:
                atom_peak_presence_shifted[atom] += 1
        
        scores_unshifted = []
        scores_shifted = []
        for i in range(0, self.molMol.GetNumAtoms()):
            if peak_presence_only:
                scores_unshifted.append(atom_peak_presence_unshifted[i])
                scores_shifted.append(atom_peak_presence_shifted[i])
            else:
                scores_unshifted.append(len(self.appearance_unshifted[i]))
                scores_shifted.append(len(self.appearance_shifted[i]))
        return scores_unshifted, scores_shifted
    
    def distance_score(self, scores_unshifted, scores_shifted, combine = True):
        scores = {}
        for i in range(0, self.molMol.GetNumAtoms()):
            if combine:
                scores[i] = scores_shifted[i] - scores_unshifted[i]
            else:
                scores[i] = scores_shifted[i]
        
        ## add min of scores to all in scores
        minScore = min(scores.values())
        for i in range(0, self.molMol.GetNumAtoms()):
            scores[i] += max(0, -minScore)
        return scores
    
    def accuracy_score(self, modificationSiteIdx, peak_presence_only = False, combine = True, getClosestMaxAtomDistance = False, getCountMax = False, getIsMax = False):
        
        scores_unshifted, scores_shifted = self.calculate_score(peak_presence_only)
        scores = self.distance_score( scores_unshifted, scores_shifted, combine)

        Max = max(scores.values())
        if Max == 0:
            return 0
        
        closestMaxAtom = 0
        distances = 0
        count = 0
        for i in range(0, self.molMol.GetNumAtoms()):
            if scores[i] == Max:
                count += 1
                distances += (self.distances[modificationSiteIdx][i]/np.amax(self.distances))
                if self.distances[modificationSiteIdx][i] < self.distances[modificationSiteIdx][closestMaxAtom]:
                    closestMaxAtom = i
        
        if getCountMax or getIsMax or getClosestMaxAtomDistance:
            res = {'score': 1 - (distances/count)}
            if getCountMax:
                res['count'] = count
            if getIsMax:
                res['isMax'] = (scores[modificationSiteIdx] == Max)
            if getClosestMaxAtomDistance:
                res['closestMaxAtomDistance'] = closestMaxAtom
        else:
            res = 1 - (distances/count)
        
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



