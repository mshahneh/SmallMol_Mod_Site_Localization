from modifinder.engines.Abtracts import AnnotationEngine
from modifinder.classes.Compound import Compound
from typing import List, Tuple
from modifinder.engines.annotation.magma import (
    fragmentation_py,
    rdkit_engine as rdkit_engine,
)
from modifinder.utilities.general_utils import mims
import networkx as nx
from rdkit import Chem


class MAGMaAnnotationEngine(AnnotationEngine):
    def __init__(self, breaks = 2, max_water_losses = 2, ionisation_mode = 0, mz_tolerance = 0.1, ppm = 40, **kwargs):
        """
        Initializes the MAGMa annotation engine

        Parameters:
            :kwargs: additional arguments
        """
        self.breaks = breaks
        self.max_water_losses = max_water_losses
        self.ionisation_mode = ionisation_mode
        self.mz_tolerance = mz_tolerance
        self.ppm = ppm
        self.args = kwargs


    def annotate(self, network: nx.DiGraph, annotate_all: bool = False, **kwargs):
        """
        Annotates the network using the MAGMa model

        Parameters:
            :network (Network): the network to be annotated
            :annotate_all (bool): whether to annotate all compounds or only the ones that have not been annotated yet
            :kwargs: additional arguments
        """
        for node in network.nodes:
            compound = network.nodes[node]["compound"]
            if compound is not None and compound.is_known:
                if annotate_all or compound.peak_fragments_map is None:
                    self.annotate_single(compound, modify_compound=True, **kwargs)


    def annotate_single(
        self, compound: Compound, modify_compound: bool = True, **kwargs
    ) -> List[Tuple[int, List[str]]]:
        """
        provides annotation for the peaks in a single compound

        Parameters:
            :compound (Compound): the compound to be annotated
            :modify_compound (bool): whether to modify the passed compound with the annotations
            :kwargs: additional arguments

        Returns:
            :Mapping[int, List[str]]: a dictionary with the indices of the peaks as keys and the list of annotations as values
        """

        if "breaks" not in kwargs:
            kwargs["breaks"] = self.breaks
        if "mz_tolerance" not in kwargs:
            kwargs["mz_tolerance"] = self.mz_tolerance
        if "ppm" not in kwargs:
            kwargs["ppm"] = self.ppm

        for item in ["breaks", "mz_tolerance", "ppm"]:
            if item not in kwargs or kwargs[item] is None:
                raise ValueError(
                    f"Missing parameters for MAGMa annotation engine {item}"
                )

        fragmentation_instance = fragmentation_py.FragmentEngine(
            compound.structure, kwargs["breaks"], max_water_losses=self.max_water_losses, ionisation_mode=self.ionisation_mode, skip_fragmentation=False, molcharge=compound.spectrum.precursor_charge
        )
        fragmentation_instance.generate_fragments()

        # generate peak to fragment map
        base_precision = 1 + kwargs["ppm"] / 1000000
        peak_fragments_map = [set() for i in range(len(compound.spectrum.mz))]
        for i in range(len(compound.spectrum.mz)):
            search_weight = compound.spectrum.mz[i] - compound.adduct_mass
            annotations = fragmentation_instance.find_fragments(
                search_weight, 0.1, base_precision, kwargs["mz_tolerance"]
            )
            for annotation in annotations:
                peak_fragments_map[i].add(annotation[0])

        if modify_compound:
            compound.peak_fragments_map = peak_fragments_map

        return peak_fragments_map


    def get_fragment_info(self, Compound: Compound, fragment):
        """
        converts a fragment to a SMILES string

        Parameters
        ----------
            Compound (Compound): the compound
            fragment: the fragment representation

        Returns
        -------
            tuple(atomlist -> List[int], edge_list -> List[Tuple[int, int]], formula -> str, smiles -> str): a tuple containing the atom list, the edge list, the formula and the SMILES string
        """
        
        atomstring = ""
        atomlist = []
        edgeList = []
        elements = dict([(e, 0) for e in mims.keys()])
        natoms = Compound.structure.GetNumAtoms()
        for atom in range(natoms):
            if 1 << atom & fragment:
                atomstring += "," + str(atom)
                atomlist.append(atom)
                elements[rdkit_engine.GetAtomSymbol(Compound.structure, atom)] += 1
                elements["H"] += rdkit_engine.GetAtomHs(Compound.structure, atom)

        for bond in Compound.structure.GetBonds():
            if bond.GetBeginAtomIdx() in atomlist and bond.GetEndAtomIdx() in atomlist:
                edgeList.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))

        formula = ""
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
            rdkit_engine.FragmentToSmiles(Compound.structure, atomlist),
        )
