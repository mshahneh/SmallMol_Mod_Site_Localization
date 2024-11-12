import modifinder.utilities.general_utils as general_utils
from modifinder.utilities.mol_utils import _get_molecule
from modifinder.utilities.spectra_utils import get_spectrum
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
from typing import Dict
from rdkit import Chem


# import modifinder as mf
from modifinder.classes.Spectrum import Spectrum
from modifinder import convert

class Compound:
    """ A class to represent a compound 

    The compound always has spectrum data including peaks, Adduct, Precursor_MZ, Charge.
    If the compound structure is known, it also has a structure property.

    Parameters
    ----------
    Main Attributes:
        id (str): The id of the compound
        spectrum (Spectrum): Spectrum Tuple representing mz, intensity, Precursor mass and Charge for mass spectrumetry data
    
    Known Compound Attributes:
        structure (Chem.Mol): The structure of the compound
        distances (dict): A dictionary of distances between every pair of atoms in the compound
        peak_fragments_map (dict): A dictionary mapping peaks to fragments

    Other Attributes:
        adduct_mass (float): The adduct mass, derived from the Adduct
        is_known (bool): A boolean indicating whether the compound is known, derived from the structure
        usi (str): The USI of the compound
        name (str): The name of the compound
        accession (str): The accession of the compound
        library_membership (str): The GNPS library membership of the compound
    
    Examples
    --------
    Create a compound by providing the necessary information:

    >>> compound = Compound(id="CCMSLIB00005435812", peaks=[[110.066925,38499.089844],[138.060638,412152.093750],[195.079575,6894530.000000],[195.200180,480874.812500],[196.082092,43027.628906]], precursor_mz=195.087, precursor_charge=1, adduct="[M+H]+", smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

    Alternatively, you can create a compound by providing a dictionary of data:

    >>> data = {
    ...     "id": "CCMSLIB00005435812",
    ...     "peaks": [[110.066925,38499.089844],[138.060638,412152.093750],[195.079575,6894530.000000],[195.200180,480874.812500],[196.082092,43027.628906]],
    ...     "Precursor_MZ": 195.087,
    ...     "Charge": 1,
    ...     "Adduct": "[M+H]+",
    ...     "Smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    ... }
    >>> compound = Compound(data)

    You can also create a compound by providing a usi:

    >>> usi = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00005435812"
    >>> compound = Compound(usi)

    or with an accession:

    >>> accession = "CCMSLIB00005435812"
    >>> compound = Compound(accession)
    
    """

    # def __init__(self, data = None, structure = None, id: str = None, spectrum: Spectrum = None,
    #             is_known: bool = None, name: str = None, peak_fragments_map: dict = None, distances: dict = None, **kwargs):
    def __init__(self, incoming_data=None, **kwargs):
        """Initialize the Compound class

        The compound class can be initialized in three different ways:
        1. By providing a dictionary of data to the *data* parameter that contains all the necessary information
        2. By providing a usi to the *data* parameter to retrieve the necessary information from GNPS
        3. By providing the necessary information as parameter

        Parameters:
        ----------
        incoming_data : input data (optional, default: None)
            The data to initialize the class with, can be a dictionary of data, a usi, or a compound object, if None, an empty compound object will be created

        kwargs : keyword arguments (optional, default: No attributes)
            Attributes to initialize the class with, if provided, they will override the attributes from the data

        See Also:
        ---------
        convert

        Examples:
        ---------

        """
        #     data (dict, str): A dictionary of data or a usi to retrieve data from GNPS
        #     structure (Chem.Mol, Smiles, InChi): The structure of the compound
        #     id (str): The id of the compound
        #     spectrum (Spectrum): an instance of Spectrum containing peak information [(mz, intensity)], precursor mass, charge, and adduct for mass spectrumetry data
        #     is_known (bool): A boolean indicating whether the compound is known, if set to False, annotators or other parts of the code will treat this compound as unknown, if not provided but the structure is provided, it will be set to True
        #     name (str): The name of the compound
        #     peak_fragments_map (dict): A dictionary mapping peaks to fragments
        #     distances (dict): A dictionary of distances between every pair of atoms in the compound, if not provided, it will be calculated from the structure
        #     **kwargs: Additional data
        #         - A use case is for the scenarios where the *data* parameter is provided, these arguments will be used to parse and clean the data
        # """
        
        # define the attributes of the class
        self.id = None
        self.spectrum = None
        self.structure = None
        self.is_known = None
        self.name = None
        self.peak_fragments_map = None
        self.distances = None
        self.usi = None
        self.adduct_mass = None

        if incoming_data is None and len(kwargs) == 0:
            return

        # attempt to initialize the class with the provided data
        if incoming_data is not None:
            convert.to_compound(incoming_data, use_object = self)
        
        # update the attributes with the provided keyword arguments
        self.update(**kwargs)



        # TODO: write setters for spectrum, structure to warn the user to update the dependent attributes
    

    def update(self, structure = None, id: str = None, spectrum: Spectrum = None, usi: str = None, 
               adduct_mass: float = None, is_known: bool = None, name: str = None,
               peak_fragments_map: dict = None, distances: dict = None, **kwargs):
        """Update the attributes of the class

        Parameters:
        ----------
        structure (Chem.Mol, Smiles, InChi): The structure of the compound
        id (str): The id of the compound
        spectrum (Spectrum): an instance of Spectrum containing peak information [(mz, intensity)], precursor mass, charge, and adduct for mass spectrumetry data
        usi (str): The USI of the compound
        adduct_mass (float): The adduct mass, derived from the Adduct
        is_known (bool): A boolean indicating whether the compound is known, if set to False, annotators or other parts of the code will treat this compound as unknown, if not provided but the structure is provided, it will be set to True
        name (str): The name of the compound
        peak_fragments_map (dict): A dictionary mapping peaks to fragments
        distances (dict): A dictionary of distances between every pair of atoms in the compound, if not provided, it will be calculated from the structure
        **kwargs: Additional data
            - A use case is for the scenarios where the *data* parameter is provided, these arguments will be used to parse and clean the data
        """
        # convert keys to lowercase
        lower_kwargs = {key.lower(): value for key, value in kwargs.items()}
        try:
            temp_structure = _get_molecule(structure, **lower_kwargs)
            self.structure = temp_structure if temp_structure is not None else self.structure
        except:
            pass
        self.id = id if id is not None else self.id
        try:
            spectrum = get_spectrum(spectrum, **lower_kwargs)
        except:
            spectrum = None
        self.spectrum = spectrum if spectrum is not None else self.spectrum
        self.usi = usi if usi is not None else self.usi
        self.adduct_mass = adduct_mass if adduct_mass is not None else self.adduct_mass
        self.is_known = is_known if is_known is not None else self.is_known
        self.name = name if name is not None else self.name
        if self.name is None:
            self.name = lower_kwargs.get('compound_name', None)
        self.peak_fragments_map = peak_fragments_map if peak_fragments_map is not None else self.peak_fragments_map
        self.distances = distances if distances is not None else self.distances

        # update the rest of the attributes
        for key, value in lower_kwargs.items():
            if key not in ['structure', 'id', 'spectrum', 'usi', 'adduct_mass', 'is_known', 'name', 'peak_fragments_map', 'distances']:
                if key not in self.spectrum.__dict__.keys():
                    setattr(self, key, value)
        
        self._parse_data()
    

    def _parse_data(self):
        """ Parse missing and verify the data of the class"""
        if self.is_known is None:
            self.is_known = (self.structure is not None)
        
        if self.distances is None and self.structure is not None:
            Chem.rdmolops.GetDistanceMatrix(self.structure)
        
        if self.adduct_mass is None and self.spectrum is not None:
            self.adduct_mass = general_utils.GetAdductMass(self.spectrum.adduct)
        
        if self.is_known is None and self.structure is not None:
            self.is_known = True
        
        # TODO: check for a valid compound
        # what is needed for a compound:
        # id, spectrum.peaks, spectrum.precursor_mass, spectrum.charge, spectrum.adduct
    

    def __str__(self):
        valuable_data_keys = ['id', 'name', 'usi']
        strings = [f"{key}: {self.__dict__[key]}" for key in valuable_data_keys if self.__dict__[key] is not None]
        joined = ', '.join(strings)

        result = f"Compound({joined}) with {len(self.spectrum.mz)} peaks"
        if self.structure is not None:
            result += f" and structure {Chem.MolToSmiles(self.structure)}"

        return result
    

    def get_meta_data(self):
        """ Get the meta data of the compound

        Returns:
            dict: A dictionary containing the meta data of the compound
        """
        description = {
            "num_peaks": len(self.spectrum.mz),
            "adduct": self.spectrum.adduct,
            "precursor_mz": self.spectrum.precursor_mass,
            "charge": self.spectrum.charge,
        }

        if self.name != None:
            description["name"] = self.name
        
        if self.is_known:
            description['num_aromatic_rings'] = rdMolDescriptors.CalcNumAromaticRings(self.structure)
            description['num_atoms'] = self.structure.GetNumAtoms()
            description['num_bonds'] = self.structure.GetNumBonds()
            description['num_rings'] = rdMolDescriptors.CalcNumRings(self.structure)

            # TODO: add peak annotation information
        
        return description


    # TODO: Implement the following methods
    # TODO: Implement the method
    def clean_peaks():
        """ Clean the peaks of the compound

        removes peaks that are small or not significant
        """
        pass

    # TODO: Implement the method
    def filter_fragments_by_atoms(self, atoms, peaks):
        """Filter the fragments by the atoms, remove fragments that do not contain at least one of the atoms
        Args:
            atoms: a list of atoms to filter the fragments
            peaks: a list of peaks to filter their fragments, if None, use all peaks
        """
        pass
    
    # TODO: Implement the method
    def calculate_peak_annotation_ambiguity(self, peaks=None):
        """Calculate the peak annotation ambiguity
        Args:
            peaks: a list of peaks to calculate the ambiguity for, if None, use all peaks
        Returns:
            ambiguity: the average number of fragments per annotated peaks
            ratio: the ratio of annotated peaks to all peaks
        """
        pass

    # TODO: Implement the method
    def calculate_annotation_entropy(self, peaks=None):
        """Calculate the entropy of the annotation
        Args:
            peaks: a list of peaks to calculate the entropy for, if None, use all peaks
        Returns:
            entropy: the entropy of the annotation
        """

    # TODO: Implement the method
    def apply_helper(self, helper_compound, shifted, unshifted, unshifted_mode = "union"):
        """use helper compound to update the peak_fragments_map
        Args:
            helper_compound: the helper compound
            shifted: array of pairs of the shifted peaks, each pair is (self_peak_index, helper_peak_index)
            unshifted: array pf pairs of the unshifted peaks, each pair is (self_peak_index, helper_peak_index)
            unshifted_mode: the mode to update the unshifted peaks, can be "union", "intersection", None
        """
