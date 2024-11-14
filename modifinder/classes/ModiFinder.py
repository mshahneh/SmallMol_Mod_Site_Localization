"""Base class for ModiFinder.

This class is used to create a ModiFinder object.
The object can be used to get information about unknown compounds in the network using the known compounds.
"""

from modifinder import convert as convert
from modifinder.classes.Compound import Compound
from modifinder.classes.EdgeDetail import EdgeDetail
from modifinder.exceptions import (
    ModiFinderNotImplementedError,
    ModiFinderNotSolvableError,
)

import networkx as nx


class ModiFinder:
    """
    Base class for ModiFinder.

    ModiFinder Can be used in two scenarios, between a known and an unknown compound. Or in a network
    consisting of known and unknown compounds.

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

    
    def solve(self, unknown: str, alignmentEngine: None, **kwargs):
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
            if alignmentEngine:
                alignmentEngine.align(self, **kwargs)
            # annotate
            # predict
            pass

    def update_edge(self, u, v, edgeDetail: EdgeDetail = None, **kwargs):
        """
        Update the edge between two compounds.

        The method will update the edge between two compounds.
        """
        if edgeDetail:
            self.network[u][v]["edgedetail"] = edgeDetail

        # update the key in the edgedetail if passed in kwargs
        for key in self.network[u][v]["edgedetail"].__dict__:
            if key in kwargs:
                setattr(self.network[u][v]["edgedetail"], key, kwargs[key])
