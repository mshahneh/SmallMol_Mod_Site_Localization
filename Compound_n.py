from rdkit import Chem
import fragmentation_py as fragmentation_py
import copy
import json
import utils_n as utils

important_arguments = ["peaks", "Adduct", "Precursor_MZ", "Charge"]
class Compound():
    """A class to represent a compound."""
    def __init__(self, data, structure=None, args={}):
        """Initialize the compound."""
        
        self.args = {"ppm": 1.01, "mz_tolerance": 0.1, "filter_peaks_method": "top_k", "filter_peaks_variable": 50}
        self.args.update(args)

        self.metadata = copy.deepcopy(data)
        for arg in important_arguments:
            if not arg in data and not arg.lower() in data:
                # print("debug: missing argument: " + arg + " in " + str(data))
                if arg != "peaks":
                    raise ValueError("Missing argument: " + arg)
                elif "peaks_json" in data:
                    self.peaks = json.loads(data["peaks_json"])
                else:
                    raise ValueError("Missing argument: " + arg)
            else:
                if arg.lower() in data:
                    setattr(self, arg, data[arg.lower()])
                    self.metadata.pop(arg.lower())
                else:
                    setattr(self, arg, data[arg])
                    self.metadata.pop(arg)
        
        self.Charge = int(self.Charge)
        self.Precursor_MZ = float(self.Precursor_MZ)
        self.Adduct = utils.parse_adduct(self.Adduct)
        self.peaks = utils.filter_peaks(self.peaks, self.args['filter_peaks_method'], self.args['filter_peaks_variable'], self.Precursor_MZ, self.Charge)

        if structure == None and "Smiles" in data:
            print("debug: using smiles")
            self.structure = Chem.MolFromSmiles(data["Smiles"])
        elif structure != None:
            if type(structure) == str:
                self.structure = Chem.MolFromSmiles(structure)
            else:
                self.structure = structure
        
        if self.structure != None:
            self.distances = Chem.rdmolops.GetDistanceMatrix(self.structure)
            self.fragments = fragmentation_py.FragmentEngine(Chem.MolToMolBlock(self.structure), 2, 2, 1, 0, 0)
            self.numFrag = self.fragments.generate_fragments()
            self.generate_peak_to_fragment_map()
    
    def __repr__(self):
        """Return a string representation of the instance."""
        obj = self.metadata
        obj.update({arg: getattr(self, arg) for arg in important_arguments})
        return json.dumps(obj)
    
    def generate_peak_to_fragment_map(self):
        self.peak_fragments_map = [set() for i in range(len(self.peaks))]
        for i in range(len(self.peaks)):
            annotations = self.fragments.find_fragments(self.peaks[i][0], 0.1, self.args['ppm'], self.args['mz_tolerance'])
            for annotation in annotations:
                self.peak_fragments_map[i].add(annotation[0])

    def apply_sirius(self, sirius):
        """Apply the sirius results to the compound and filter the possible annotations"""
        for i, peak in enumerate(self.peaks):
            index = utils.find_mz_in_sirius(sirius['fragments'], peak[0], self.args['mz_tolerance'])
            if index == -1:
                continue
            helper_peak_formula = sirius['fragments'][index]['molecularFormula']
            possibilites = set()
            for frag_id in self.peak_fragments_map[i]:
                molSubFormula = self.fragments.get_fragment_info(frag_id, 0)[2]
                if utils.is_submolecule(molSubFormula, helper_peak_formula) and utils.is_submolecule(helper_peak_formula, molSubFormula):
                    possibilites.add(frag_id)
            if len(possibilites) > 0:
                self.peak_fragments_map[i] = possibilites
    
    
    def apply_helper(self, helper_compound):
        # TODO: implement
        pass

    def apply_iceberg(self, iceberg):
        # TODO: implement
        pass
