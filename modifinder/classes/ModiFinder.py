import modifinder as mf
import networkx as nx

class ModiFinder():
    def __init__(self, compounds: list, unknowns: list, network: nx.Graph = None, **kwargs):
        """
        """
        self.network = None
        self.unknonws = None
        self.args = kwargs


        # define the attributes of the class
        if not isinstance(unknonws, list):
            unknonws = [unknowns]
        
        if self.network is not None:
            ## if the network is passed, then the compounds and unknowns are not needed 
            self.network = network
            if unknowns is not None:
                self.unknowns = unknowns
            else:
                self.unknowns = []
                for node in self.network.nodes:
                    if not self.network.nodes[node].get('compound').is_known:
                        self.unknowns.append(node)
        else:
            ## if the network is not passed, then the compounds and unknowns are used to create the network
            self.network = nx.Graph()
            for compound in compounds:
                if not isinstance(compound, mf.Compound):
                    if isinstance(compound, str):
                        compound = mf.Compound(usi=compound)
                    else:
                        compound = mf.Compound(data=compound)
                self.network.add_node(compound.id, compound=compound)
            self.unknowns = unknowns

    # check input is in correct format
    def _verify_object():
        # check all the compounds have ids
        pass


    
    def add_helpers():
        pass

    def add_edge():
        pass
