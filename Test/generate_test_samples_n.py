"""   
Compound_Item format:
    {
        "peaks":[],
        "Adduct": str,
        "Precursor_MZ": float,
        "Charge": int,
        "Smiles": str,
        "sirius_path": str,
    }

Expected_Result_Item format:
    {
        peak_fragment_map: [[int]],
        fragments: [],
    }

Match_Item format:
    {
        "type": str ("synthetic" or "database"),
        "library": str,
        "main": Compound_Item,
        "helpers": [Compound_Item],
        "modified": Compound_Item,
        "expected_result": Expected_Result_Item
    }
"""

from urllib.request import urlopen
from IPython.display import SVG
import matplotlib.pyplot as plt
from rdkit import Chem
import pickle
import os
import sys
import json
from tqdm import tqdm
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)

import library_downloader as library_downloader

def generate_compound_item(data, libpath, id):
    helperDirectory = os.path.join(libpath,"nf_output/fragmentationtrees/")
    result = {
            "peaks": json.loads(data["peaks_json"]),
            "Adduct": data["Adduct"],
            "Precursor_MZ": data["Precursor_MZ"],
            "Charge": data["Charge"],
            "Smiles": data["Smiles"],
        }
    try:
        molSirius = json.load(open(os.path.join(helperDirectory, id + "_fragmentationtree.json")))
        result["sirius"] = os.path.join(helperDirectory, id + "_fragmentationtree.json")
    except:
        pass
    return result


def generate_test_from_library(lib, count):
    # if directory does not exist, create it
    abs_path = os.path.abspath("../data/libraries")
    location = os.path.join(abs_path, lib)
    if not os.path.exists(location):
        url = "https://gnps-external.ucsd.edu/gnpslibrary/" + lib + ".json"
        library_downloader.download(url, location, 0.5, 0.1)
    
    # load data
    with open(os.path.join(location,"data_dict_filtered.pkl"), "rb") as f:
        data_dict_filtered = pickle.load(f)
    with open(os.path.join(location,"matches.pkl"), "rb") as f:
        matches = pickle.load(f)
    
    # compute helper pairs
    helpers = dict()
    for match in matches[1]:
        if match[0] not in helpers:
            helpers[match[0]] = []
        helpers[match[0]].append(match[1])
    
    results = []
    # generate test samples
    for match in matches[1]:
        if len(results) >= count:
            break

        if (data_dict_filtered[match[0]]["Adduct"] != "M+H" or data_dict_filtered[match[0]]["Adduct"] != data_dict_filtered[match[1]]["Adduct"]):
            continue

        main = generate_compound_item(data_dict_filtered[match[1]], location, match[1])
        helpers_compound_items = []
        for helper in helpers.get(match[0], []):
            helpers_compound_items.append(generate_compound_item(data_dict_filtered[helper], location, helper))
        modified = generate_compound_item(data_dict_filtered[match[0]], location, match[0])
        results.append({
            "type": "database", 
            "library": lib,
            "main": main,
            "helpers": helpers_compound_items,
            "modified": modified,
            "expected_result": None
        })
    
    return results
    

def generate_synthetic_test(count):
    return []

def main(libs, lib_count, synthetic_count):
    """
        Generate test samples
        output:
            [Match_Item] -> array of Match_Item
    """
    results = []
    for lib in libs:
        results.extend(generate_test_from_library(lib, lib_count))
    results.extend(generate_synthetic_test(synthetic_count))
    return results
    
if __name__ == '__main__':
    libs = ["BERKELEY-LAB"]
    results = main(libs, 5, 1)
    with open("test_samples.pkl", "wb") as f:
        pickle.dump(results, f)