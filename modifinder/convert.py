import warnings
import modifinder as mf
import modifinder.utilities.network as network
import json
import modifinder.utilities.gnps_types as gt
from copy import deepcopy

def to_compound(data, use_object=None):
    """Make a Compound object from the data
    
    Parameters
    ----------
    data: object to be converted

        Current supported types are:
         Compound object
         USI string
         dictionary-of-data
    
    use_object: object, optional
        If a Compound object is passed, this object will be used to create the new object.
    """

    # Compound Object
    if hasattr(data, "spectrum"):
        try:
            if use_object:
                compound = use_object
                compound.clear()
                compound.update(compound_to_dict(data))
            else:
                compound = mf.Compound()
                copied_data = deepcopy(compound_to_dict(data))
                compound.update(**copied_data)
            return compound

        except Exception as err:
            raise mf.ModiFinderError(f"Input data is not a valid Compound object.") from err
    
    # USI
    if isinstance(data, str):
        try:
            data = network.get_data(data)
            if use_object:
                compound = use_object
                compound.clear()
                compound.update(**data)
            else:
                compound = mf.Compound()
                compound.update(**data)
            return compound

        except Exception as err:
            raise mf.ModiFinderError(f"Input data is not a valid USI string.") from err
    
    # Dictionary
    if isinstance(data, dict):
        try:
            if use_object:
                compound = use_object
                compound.clear()
                compound.update(**data)
            else:
                compound = mf.Compound()
                compound.update(**data)
            return compound
        except Exception as err:
            raise mf.ModiFinderError(f"Input data is not a valid dictionary.") from err
        

def compound_to_dict(compound):
    """Convert a Compound object to a dictionary"""
    return compound.__dict__



def convert_to_universal_key(key: str) -> str:
    """
    Convert different types of keys to universal keys.
    This function standardizes various key names to a universal format. 

    Args:
        :key (str): The key to be converted.
    
    Returns:
        :str: The converted key.
    """
    key = key.lower()
    key = key.replace(" ", "_")
    return gt.gnps_keys_mapping.get(key, key)
    
# TODO: use machine learning prepared data instead of hardcoding
def parse_data_to_universal(data):
    """
    Parse the data to a universal format.

    This function takes a dictionary of data and converts it into a universal format.
    It processes specific keys like "peaks_json" and "Charge" differently, and attempts
    to convert other values to floats. If the conversion to float is successful and the
    key is "Charge", it further converts the value to an integer.

    Args:
        :data (dict): The input data dictionary to be parsed.

    Returns:
        :dict: A dictionary with keys converted to a universal format and values processed
              accordingly.
    """

    res = {}
    for key, value in data.items():
        converted_key = convert_to_universal_key(key)
        if key == "peaks_json":
            res['peaks'] = json.loads(value)
        elif converted_key == "adduct":
            res[converted_key] = gt.adduct_mapping.get(value, value)
        else:
            try:
                value = float(value)
                if key in ["charge", "ms_level"]:
                    value = int(value)
            except:
                pass
            res[converted_key] = value
    return res