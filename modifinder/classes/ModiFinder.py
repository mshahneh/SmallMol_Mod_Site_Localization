"""Base class for ModiFinder.

This class is used to create a ModiFinder object.
The object can be used to get information about unknown compounds in the network using the known compounds.
"""

from modifinder import convert as convert
from modifinder.classes.Compound import Compound
from modifinder.classes.EdgeDetail import EdgeDetail, MatchType
from modifinder.engines.Abtracts import AlignmentEngine, AnnotationEngine
from modifinder.engines.annotation.MAGMaAnnotationEngine import MAGMaAnnotationEngine
from modifinder.engines.alignment.CosineAlignmentEngine import CosineAlignmentEngine

from modifinder.exceptions import (
    ModiFinderNotImplementedError,
    ModiFinderNotSolvableError,
)

import networkx as nx
import numpy as np
from modifinder.utilities import general_utils as general_utils


class ModiFinder:
    """
    Base class for compound network.

    ModiFinder Can be used in two scenarios
    1. between a known and an unknown compound. 
    2. in a network consisting of known and unknown compounds.

    In both scenarios, the ModiFinder object stores a directed graph where the nodes are compounds (known and unknown)
    and edges are directed (from smaller compound to the larger compound) holding the relationships between the pair (alignment, similarity, etc). The object also stores
    the unknown compounds in the network. The object can be used to get information about the unknown compound.

    Parameters
    ----------
    network : nx.DiGraph
        A networkx graph object with the nodes are identified by the compound ids and the "compound"
        attibute of a node is a Compound object. The edges are the relationships between the compounds
        and the "edgedetail" attribute of an edge is an EdgeDetail object.

    unknowns : list
        A list of compound ids that are unknown in the network.

    See Also
    --------
    Compound
    EdgeDetail
    Engines

    Examples
    --------
    """

    def __init__(
        self,
        knownCompond: Compound = None,
        unknownCompound: Compound = None,
        edgeDetail: EdgeDetail = None,
        network: nx.DiGraph = None,
        networkUnknowns: list = None,
        should_align: bool = True,
        alignmentEngine: AlignmentEngine = None,
        should_annotate: bool = True,
        annotationEngine: AnnotationEngine = None,
        **kwargs,
    ):
        """
        Initialize the ModiFinder object.

        **Use Case  1: a pair of known and unknown compound**

        **Use Case  2: a network of known and unknown compounds**

        Parameters
        ----------
        knownCompond : known compound (not optional in Use Case 1, ignored in Use Case 2)
            Data to create a known compound. The data can be a Compound object, or any data that can be converted to a Compound object.

        unknownCompound : unknown compound (not optional in Use Case 1, ignored in Use Case 2)
            Data to create an unknown compound. The data can be a Compound object, or any data that can be converted to a Compound object.

        edgeDetail : EdgeDetail object from smaller to larger compound (optional in Use Case 1, ignored in Use Case 2)
            the orientation of the match must be from the smaller compound to the larger compound.

        network : nx.Graph (if passed, then Use Case 2 is used, if None, then Use Case 1 is used)
            A networkx graph object with the nodes are identified by the compound ids and the "compound"
            attibute of a node is a *Compound* object. The edges are the relationships between the compounds
            and the "edgedetail" attribute of an edge is an *EdgeDetail* object.

        networkUnknowns : list (ignored in Use Case 1, optional in Use Case 2)
            A list of compound ids that are unknown in the network, if not passed, the unknowns will be driven from 'is_known' attribute
            of the compounds in the network.


        See Also
        --------
        Compound
        EdgeDetail
        Engines

        Examples
        --------

        """

        self.network = None
        self.unknonws = None
        self.args = kwargs

        if network:
            self.network = network

            if networkUnknowns:
                self.unknowns = networkUnknowns
            else:
                self.unknowns = [
                    node
                    for node in network.nodes
                    if not network.nodes[node]["compound"].is_known
                ]

        else:
            self.network = nx.DiGraph()
            knownCompond = convert.to_compound(data=knownCompond)
            unknownCompound = convert.to_compound(data=unknownCompound)
            self.network.add_node(knownCompond.id, compound=knownCompond)
            self.network.add_node(unknownCompound.id, compound=unknownCompound)

            if (
                knownCompond.spectrum.precursor_mz
                < unknownCompound.spectrum.precursor_mz
            ):
                self.network.add_edge(
                    knownCompond.id, unknownCompound.id, edgedetail=edgeDetail
                )
            else:
                self.network.add_edge(
                    unknownCompound.id, knownCompond.id, edgedetail=edgeDetail
                )

            self.unknowns = [unknownCompound.id]
        
        if should_align:
            if alignmentEngine is None:
                alignmentEngine = CosineAlignmentEngine()
            alignmentEngine.align(self.network, **kwargs)
        
        if should_annotate:
            if annotationEngine is None:
                annotationEngine = MAGMaAnnotationEngine()
            annotationEngine.annotate(self.network, annotate_all = True, **kwargs)
    
    
    def align():
        pass
    
    def annotate():
        pass
    
    def align_and_annotate():
        pass

    
    def solve(self, unknown: str, **kwargs):
        """
        Solve the network.

        The method will solve the network for an unknown compound
        
        Parameters
        ----------
        unknown : str
            The id of the unknown compound in the network.
        
        alignmentEngine : AlignmentEngine
            The alignment engine to use to align the unknown compound with the known compounds in the network.
        """

        if self.network.in_degree(unknown) + self.network.out_degree(unknown) == 0:
            raise ModiFinderNotSolvableError(
                f"{unknown} is not solveable due to lack of known compounds"
            )

        if self.network.in_degree(unknown) + self.network.out_degree(unknown) > 1:
            # TODO: implement the logic to solve the network when there are multiple known compounds
            raise ModiFinderNotImplementedError(
                f"{unknown} is not solveable due to multiple known compounds"
            )

        if self.network.in_degree(unknown) == 1:
            # TODO: implement the logic to solve the network when there is only one known compound
            pass


    def update_edge(self, u, v, edgeDetail: EdgeDetail = None, **kwargs):
        """
        Update the edge between two compounds.

        The method will update the edge between two compounds.
        """
        if edgeDetail:
            self.network[u][v]["edgedetail"] = edgeDetail
        else:
            if not u.spectrum.precursor_mz <= v.spectrum.precursor_mz:
                raise ValueError(
                    f"the edge between {u.id} and {v.id} must be from the smaller compound to the larger compound"
                )
            else:
                self.network[u][v]["edgedetail"] = EdgeDetail()
                

        # update the key in the edgedetail if passed in kwargs
        for key in self.network[u][v]["edgedetail"].__dict__:
            if key in kwargs:
                setattr(self.network[u][v]["edgedetail"], key, kwargs[key])
    
    
    def generate_probabilities(self, known_id = None, unknown_id = None, shifted_only = False, CI = False, CPA = True, CFA = True, CPE = True):
        
        if unknown_id is None:
            if len(self.unknowns) > 1:
                raise ValueError("More than one unknown compound found in the network. Please specify the unknown compound id")
            unknown_id = self.unknowns[0]
        
        if known_id is None and unknown_id is not None:
            neighbors = list(self.network.predecessors(unknown_id)) + list(self.network.successors(unknown_id))
            if len(neighbors) > 1:
                raise ValueError("More than one known compound found in the network. Please specify the known compound id")
            known_id = neighbors[0]
        
        if known_id is None or unknown_id is None:
            raise ValueError("Both known and unknown compound ids must be specified")
        
        # check if unknown is connected to known
        if self.network.has_edge(known_id, unknown_id):
            edgeDetail = self.network[known_id][unknown_id]["edgedetail"]
            if edgeDetail is None:
                raise ValueError(f"Edge between {known_id} and {unknown_id} does not have edge details")
        
            shifted_peaks_in_known = [match.first_peak_index for match in edgeDetail.matches if match.match_type == MatchType.shifted]
            unshifted_peaks_in_known = [match.first_peak_index for match in edgeDetail.matches if match.match_type == MatchType.unshifted]
        elif self.network.has_edge(unknown_id, known_id):
            edgeDetail = self.network[unknown_id][known_id]["edgedetail"]
            if edgeDetail is None:
                raise ValueError("No edge detail found between the known and unknown compounds")    
            shifted_peaks_in_known = [match.second_peak_index for match in edgeDetail.matches if match.match_type == MatchType.shifted]
            unshifted_peaks_in_known = [match.second_peak_index for match in edgeDetail.matches if match.match_type == MatchType.unshifted]
        else:
            raise ValueError("No edge found between the known and unknown compounds")
        known_compound = self.network.nodes[known_id]["compound"]
        positive_contributions = known_compound.calculate_contributions(shifted_peaks_in_known, CI=CI, CPA=CPA, CFA=CFA, CPE=CPE)
        if not shifted_only:
            negative_contributions = known_compound.calculate_contributions(unshifted_peaks_in_known, CI=CI, CPA=CPA, CFA=CFA, CPE=CPE)
        else:
            negative_contributions = [0 for i in range(len(known_compound.structure.GetAtoms()))]
        
        probabilities = np.zeros(len(known_compound.structure.GetAtoms()))
        for i in range(len(positive_contributions)):
            probabilities[i] = positive_contributions[i] - negative_contributions[i]
        
        if np.min(probabilities) < 0:
            probabilities = probabilities - np.min(probabilities)
            # mask only the atoms with positive contribution
            # probabilities = [probabilities[i] if positive_contributions[i] > 0 else 0 for i in range(len(probabilities))]
        if np.sum(probabilities) > 0:
            probabilities = probabilities / np.sum(probabilities)

        probabilities = general_utils.power_prob(probabilities)

        return probabilities
