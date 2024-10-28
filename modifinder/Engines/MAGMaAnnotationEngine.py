from modifinder.Engines.Abtracts import AnnotationEngine
from modifinder.Compound import Compound
from typing import List, Tuple
from modifinder.Engines.magma import fragmentation_py
from modifinder.utilities.general_utils import mims
from rdkit import Chem

def fragment_to_smiles(mol, atomlist):
    emol = Chem.EditableMol(mol)
    for atom in reversed(range(mol.GetNumAtoms())):
        if atom not in atomlist:
            emol.RemoveAtom(atom)

    frag = emol.GetMol()
    return Chem.MolToSmiles(frag)


class MAGMaAnnotationEngine(AnnotationEngine):
    def __init__(self, **kwargs):
        """
        Initializes the MAGMa annotation engine

        Parameters:
            :kwargs: additional arguments
        """
        self.args = kwargs
        pass


    def annotate(self, network, **kwargs):
        """
        Annotates the network using the MAGMa model

        Parameters:
            :network (Network): the network to be annotated
            :kwargs: additional arguments
        """
        pass


    def annotate_single(self, compound: Compound,  modify_compound: bool = True, **kwargs) -> List[Tuple[int, List[str]]]:
        """
        provides annotation for the peaks in a single compound

        Parameters:
            :compound (Compound): the compound to be annotated
            :modify_compound (bool): whether to modify the passed compound with the annotations
            :kwargs: additional arguments
        
        Returns:
            :Mapping[int, List[str]]: a dictionary with the indices of the peaks as keys and the list of annotations as values
        """

        if 'breaks' not in kwargs:
            kwargs['breaks'] = self.args["breaks"]
        if 'mz_tolerance' not in kwargs:
            kwargs['mz_tolerance'] = self.args["mz_tolerance"]
        if 'ppm' not in kwargs:
            kwargs['ppm'] = self.args["ppm"]
        
        for item in ['breaks', 'mz_tolerance', 'ppm']:
            if item not in kwargs or kwargs[item] is None:
                raise ValueError(f"Missing parameters for MAGMa annotation engine {item}")
        
        fragmentation_instance = fragmentation_py.FragmentEngine(compound.structure, kwargs['breaks'], 2, 0, 0, 0)
        fragmentation_instance.generate_fragments()

        # generate peak to fragment map
        base_precision = 1 + kwargs["ppm"] / 1000000
        peak_fragments_map = [set() for i in range(len(compound.peaks))]
        for i in range(len(compound.peaks)):
            search_weight = compound.peaks[i][0] - compound.Adduct_Mass
            annotations = fragmentation_instance.find_fragments(search_weight, 0.1, base_precision, kwargs["mz_tolerance"])
            for annotation in annotations:
                peak_fragments_map[i].add(annotation[0])
        
        if modify_compound:
            compound.peak_fragments_map = peak_fragments_map
        
        return peak_fragments_map
    

    def get_fragment_info(self, Compound: Compound, fragment):
        atomstring = ''
        atomlist = []
        edgeList = []
        elements = dict([(e,0) for e in mims.keys()])
        natoms = Compound.structure.GetNumAtoms()
        for atom in range(natoms):
            if 1 << atom & fragment:
                atomstring += ',' + str(atom)
                atomlist.append(atom)
                elements[Chem.GetAtomSymbol(Compound.structure, atom)] += 1
                elements['H'] += Chem.GetAtomHs(Compound.structure, atom)
        
        for bond in self.bonds:
            if bond[0] in atomlist and bond[1] in atomlist:
                edgeList.append((bond[0], bond[1]))

        formula = ''
        for el in mims.keys():
            nel = elements[el]
            if nel > 0:
                formula += el
            if nel > 1:
                formula += str(nel)

        return (
            atomlist,
            edgeList,
            formula,
            fragment_to_smiles(Compound.structure, atomlist))
    
