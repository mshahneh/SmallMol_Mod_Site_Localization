import collections
import json

SpectrumTuple = collections.namedtuple(
    "SpectrumTuple", ["precursor_mz", "precursor_charge", "mz", "intensity"]
)


def convert_to_SpectrumTuple(peaks: list, precursor_mz: float, precursor_charge: int) -> SpectrumTuple:
    """
    Converts the peaks to SpectrumTuple.
    :param peaks: List of tuples (mz, intensity) or List of List [mz, intensity]
    :param precursor_mz: Precursor m/z
    :param precursor_charge: Precursor charge
    """
    if not peaks:
        return None
    if len(peaks) > 0 and (not isinstance(peaks[0], tuple) and not isinstance(peaks[0], list)):
        raise ValueError("Peaks should be a list of tuples (mz, intensity)")
    
    res = {}
    res['precursor_charge'] = precursor_charge
    res['precursor_mz'] = precursor_mz
    res['mz'] = []
    res['intensity'] = []
    for peak in peaks:
        res['mz'].append(peak[0])
        res['intensity'].append(peak[1])
    
    return SpectrumTuple(**res)

def convert_to_SpectrumTuple_seprated(mz: list, intensity: list, precursor_mz: float, precursor_charge: int) -> SpectrumTuple:
    """
    Converts the peaks to SpectrumTuple.
    :param mz: List of mz values
    :param intensity: List of intensity values
    :param precursor_mz: Precursor m/z
    :param precursor_charge: Precursor charge
    """
    if not mz:
        return None
    if len(mz) != len(intensity):
        raise ValueError("Length of mz and intensity should be same")
    
    res = {}
    res['precursor_charge'] = precursor_charge
    res['precursor_mz'] = precursor_mz
    res['mz'] = mz
    res['intensity'] = intensity
    
    return SpectrumTuple(**res)

def convert_to_universal_key(key: str) -> str:
    """
    Convert different types of keys to universal keys
    :param key: str
    """
    if key == "precursor_mz":
        return "Precursor_MZ"
    if key == "smiles" or key == "SMILES":
        return "Smiles"
    if key == "charge":
        return "Charge"
    if key == "adduct":
        return "Adduct"
    # TODO: Add more keys
    return key
    
def parse_data_to_universal(data):
    """
    Parse the data to universal format
    :param data: dict
    """
    res = {}
    for key, value in data.items():
        if key == "peaks_json":
            res['peaks'] = json.loads(value)
        else:
            try:
                value = float(value)
                if key == "Charge":
                    value = int(value)
            except:
                pass
            res[convert_to_universal_key(key)] = value
    return res