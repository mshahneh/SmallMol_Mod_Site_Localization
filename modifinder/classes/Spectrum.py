from typing import Dict
import json
from modifinder.utilities.gnps_types import adduct_mapping

class Spectrum:
    """A class to represent a spectrum.
    Parameters:
    ----------
    mz: list
        A list of m/z values.
    intensity: list
        A list of intensity values.
    precursor_mz: float
        The precursor m/z value.
    precursor_charge: int
        The precursor charge.
    adduct: str
        The adduct.
    ms_level: int
        The ms level, default is 2.
    instrument: str, optional
        The instrument used.
    ms_mass_analyzer: str, optional
        The mass analyzer used.
    ms_dissociation_method: str, optional
        The dissociation method used.
    spectrum_id: str, optional
    """
    def __init__(self, incoming_data:Dict=None, **kwargs):
        """Constructor for the Spectrum class.

        the spectrum class can be initialized with a dictionary of data or with the individual values.

        Args:
            incoming_data (dict): A dictionary of data to initialize the Spectrum object.
        """
        self.mz = None
        self.intensity = None
        self.precursor_mz = None
        self.precursor_charge = None
        self.adduct = None
        self.ms_level = None
        self.instrument = None
        self.ms_mass_analyzer = None
        self.ms_dissociation_method = None
        self.spectrum_id = None

        if incoming_data is not None:
            self.update(**incoming_data)

        self.update(**kwargs)


    def update(self, peaks = None, peaks_json = None, mz=None, intensity=None, precursor_mz=None, precursor_charge=None, 
               adduct=None, ms_level=None, instrument=None, ms_mass_analyzer=None, 
               ms_dissociation_method=None, spectrum_id=None, **kwargs):
        """Update the Spectrum object with the given values.

        Args:
            mz (list): A list of m/z values.
            intensity (list): A list of intensity values.
            precursor_mz (float): The precursor m/z value.
            precursor_charge (int): The precursor charge.
            adduct (str): The adduct.
            ms_level (int): The ms level.
            instrument (str): The instrument used.
            ms_mass_analyzer (str): The mass analyzer used.
            ms_dissociation_method (str): The dissociation method used.
            spectrum_id (str): The spectrum id.
        """
        if peaks_json is not None:
            peaks = json.loads(peaks_json)
        if peaks is not None:
            self.mz = [peak[0] for peak in peaks]
            self.intensity = [peak[1] for peak in peaks]
        self.mz = mz if mz is not None else self.mz
        self.intensity = intensity if intensity is not None else self.intensity
        self.precursor_mz = precursor_mz if precursor_mz is not None else self.precursor_mz
        self.precursor_charge = precursor_charge if precursor_charge is not None else self.precursor_charge
        self.adduct = adduct_mapping[adduct] if adduct is not None else self.adduct
        self.ms_level = ms_level if ms_level is not None else self.ms_level
        self.instrument = instrument if instrument is not None else self.instrument
        self.ms_mass_analyzer = ms_mass_analyzer if ms_mass_analyzer is not None else self.ms_mass_analyzer
        self.ms_dissociation_method = ms_dissociation_method if ms_dissociation_method is not None else self.ms_dissociation_method
        self.spectrum_id = spectrum_id if spectrum_id is not None else self.spectrum_id


    def __str__(self):
        object_dict = self.__dict__
        to_delete = [keys for keys in object_dict.keys() if object_dict[keys] is None]
        for key in to_delete:
            del object_dict[key]
        return json.dumps(object_dict, indent=4)

