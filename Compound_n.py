from rdkit import Chem
from msbuddy import assign_subformula
from . import fragmentation_py as fragmentation_py
import copy
import json
from . import utils_n as utils
from . import handle_network as handle_network
import re

important_arguments = ["peaks", "Adduct", "Precursor_MZ", "Charge"]


class Compound:
    """A class to represent a compound."""

    def __init__(self, data, structure=None, args={}):
        """Initialize the compound."""

        self.args = {
            "ppm": 40,
            "mz_tolerance": 0.1,
            "filter_peaks_method": "intensity",
            "filter_peaks_variable": 0.01,
        }
        self.args.update(args)
        for arg in list(args.keys()):
            # if its a string and can be converted to float, convert it
            if type(args[arg]) == str:
                try:
                    self.args[arg] = float(args[arg])
                except:
                    self.args[arg] = args[arg]

        if type(data) == str:
            data = handle_network.getDataFromUsi(data)

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
        self.peaks = utils.filter_peaks(
            self.peaks,
            self.args["filter_peaks_method"],
            self.args["filter_peaks_variable"],
            self.Precursor_MZ,
            self.Charge,
        )

        if structure == None and "Smiles" in data:
            self.structure = Chem.MolFromSmiles(data["Smiles"])
        elif structure != None:
            if type(structure) == str:
                if len(structure) > 0:
                    self.structure = Chem.MolFromSmiles(structure)
                else:
                    self.structure = None
            else:
                self.structure = structure
        else:
            self.structure = None

        if self.structure != None:
            self.distances = Chem.rdmolops.GetDistanceMatrix(self.structure)
            breaks = 5
            if (self.structure.GetNumAtoms() > 30):
                breaks = 4
            if (self.structure.GetNumAtoms() > 50):
                breaks = 3
            if (self.structure.GetNumAtoms() > 80):
                breaks = 2
            self.fragments = fragmentation_py.FragmentEngine(
                Chem.MolToMolBlock(self.structure), breaks, 2, 1, 0, 0
            )
            self.numFrag = self.fragments.generate_fragments()
            self.generate_peak_to_fragment_map()

    def __repr__(self):
        """Return a string representation of the instance."""
        obj = self.metadata
        obj.update({arg: getattr(self, arg) for arg in important_arguments})
        return json.dumps(obj)

    def generate_peak_to_fragment_map(self):
        base_precision = 1 + self.args["ppm"] / 1000000
        self.peak_fragments_map = [set() for i in range(len(self.peaks))]
        for i in range(len(self.peaks)):
            annotations = self.fragments.find_fragments(
                self.peaks[i][0], 0.1, base_precision, self.args["mz_tolerance"]
            )
            for annotation in annotations:
                self.peak_fragments_map[i].add(annotation[0])

    def apply_sirius(self, sirius):
        """Apply the sirius results to the compound and filter the possible annotations"""
        for i, peak in enumerate(self.peaks):
            index = utils.find_mz_in_sirius(
                sirius["fragments"], peak[0], self.args["mz_tolerance"], self.args["ppm"]
            )
            if index == -1:
                continue
            helper_peak_formula = sirius["fragments"][index]["molecularFormula"]
            possibilites = set()
            for frag_id in self.peak_fragments_map[i]:
                molSubFormula = self.fragments.get_fragment_info(frag_id, 0)[2]
                if utils.is_submolecule(
                    molSubFormula, helper_peak_formula
                ) and utils.is_submolecule(helper_peak_formula, molSubFormula):
                    possibilites.add(frag_id)

            self.peak_fragments_map[i] = possibilites

    def apply_filter(self, filter_method, filter_value):
        self.peaks = utils.filter_peaks(
            self.peaks, filter_method, filter_value, self.Precursor_MZ, self.Charge
        )

    def remove_large_peaks(self, weight=None):
        if weight == None:
            weight = self.Precursor_MZ - 1.6

        peaks = []
        for i in range(len(self.peaks)):
            if self.peaks[i][0] < weight:
                peaks.append(self.peaks[i])
        self.peaks = peaks

    def get_peak_index(self, peak_weight, args={}):
        args.update(self.args)
        ind = []
        for i in range(len(self.peaks)):
            diff = abs(self.peaks[i][0] - peak_weight)
            if diff < args["mz_tolerance"] and diff / peak_weight * 1000000 < args["ppm"]:
                ind.append(i)
        return ind


    def filter_fragments_by_atoms(self, atoms):
        """Filter the fragments by the atoms, remove fragments that do not contain the atoms"""
        for i in range(len(self.peak_fragments_map)):
            updated_fragments = set()
            for fragment in self.peak_fragments_map[i]:
                for atom in atoms:
                    if 1 << atom & fragment:
                        updated_fragments.add(fragment)
                        break
            self.peak_fragments_map[i] = updated_fragments


    def vote_for_fragments(self, peak_index, fragments, voter="none", extra_args={}):
        """Vote for the fragments, update the peak_fragments_map with the intersection of the fragments, and extend the intersection to other peaks if extend_to_other_peaks is true"""
        if peak_index == -1:
            print("in vote for fragment, Peak not found")
            return

        # get the intersection of the fragments
        intersection = self.peak_fragments_map[peak_index].intersection(fragments)
        self.peak_fragments_map[peak_index] = intersection
        extra_args.update(self.args)
        # if kwargs has extend_to_other_peaks, extend the intersection to other peaks
        if "extend_to_other_peaks" in extra_args and extra_args["extend_to_other_peaks"]:
            # find atoms in the intersection
            intersect_atoms = set()
            for atom in range(len(self.structure.GetAtoms())):
                exists = True
                for fragment in intersection:
                    if not 1 << atom & fragment:
                        exists = False
                        break
                if exists:
                    intersect_atoms.add(atom)
            
            self.filter_fragments_by_atoms(intersect_atoms)
        return
    

    def calculate_peak_annotation_ambiguity(self, peaks=None):
        """Calculate the peak annotation ambiguity
        Args:
            peaks: a list of peaks to calculate the ambiguity for, if None, use all peaks
        Returns:
            ambiguity: the average number of fragments per annotated peaks
            ratio: the ratio of annotated peaks to all peaks
        """
        if peaks == None:
            peaks = [i for i in range(len(self.peaks))]
        ambiguity = 0
        annotated_peaks = 0
        for peak in peaks:
            if len(self.peak_fragments_map[peak]) > 0:
                annotated_peaks += 1
                ambiguity += len(self.peak_fragments_map[peak])
        if annotated_peaks == 0:
            return -1, 0
        if len(peaks) == 0:
            return -1, -1
        return ambiguity / annotated_peaks, annotated_peaks / len(peaks)
    
    def apply_msbuddy(self):
        if not "M+H" in self.Adduct:
            raise ValueError("Adduct not supported")
        main_compound_formula = Chem.rdMolDescriptors.CalcMolFormula(self.structure)
        peak_mz = [peak[0] for peak in self.peaks]
        subformla_list = assign_subformula(peak_mz,
                                        precursor_formula=main_compound_formula, adduct="[M+H]+",
                                        ms2_tol=self.args["ppm"], ppm=True, dbe_cutoff=-1.0)
        
        for i in range(len(self.peaks)):
            possibilites = set()
            for subformula in subformla_list[i].subform_list:
                formula = subformula.formula
                
                # remove one H from the formula
                pattern = r'([A-Z][a-z]*)(\d*)'
                matches = re.findall(pattern, formula)  # Find all matches in the formula
                formula = ""
                for match in matches:
                    if match[0] == "H":
                        count = match[1]
                        if count:
                            count = int(count)
                        else:
                            count = 1
                        count -= 1
                        if count > 1:
                            formula += "H" + str(count)
                        elif count == 1:
                            formula += "H"
                    else:
                        formula += match[0]
                        if match[1]:
                            formula += match[1]
                
                # find the fragments that contains the formula
                for frag_id in self.peak_fragments_map[i]:
                    molSubFormula = self.fragments.get_fragment_info(frag_id, 0)[2]
                    if utils.is_submolecule(molSubFormula, formula) and utils.is_submolecule(formula, molSubFormula):
                        possibilites.add(frag_id)
            self.peak_fragments_map[i] = possibilites

    def apply_helper(self, helper_compound):
        # TODO: implement
        pass

    def apply_iceberg(self, iceberg):
        # TODO: implement
        pass
