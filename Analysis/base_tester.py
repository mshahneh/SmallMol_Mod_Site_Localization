from urllib.request import urlopen
from IPython.display import SVG
import matplotlib.pyplot as plt
from rdkit import Chem
from tqdm import tqdm
import pandas as pd
import xlsxwriter
import argparse
import pickle
import numpy as np
import json
import os
import sys
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)

import library_downloader as library_downloader
import multiprocessing as mp
from functools import partial


from Compound_n import Compound
from ModificationSiteLocator import ModificationSiteLocator


def load(library, accepted_adduct = None, count = None):
    print("Loading library: " + library)
    if not os.path.exists( os.path.join("../data/libraries", library)):
        url = "https://gnps-external.ucsd.edu/gnpslibrary/" + library + ".json"
        location = "../data/libraries/" + library + "/"
        library_downloader.download(url, location, 0.5, 0.1)

    with open(os.path.join("../data/libraries", library, "data_dict_filtered.pkl"), "rb") as f:
        data_dict_filtered = pickle.load(f)

    # load matches
    with open(os.path.join("../data/libraries", library, "matches.pkl"), "rb") as f:
        matches = pickle.load(f)

    # load cachedStructures_filtered
    with open(os.path.join("../data/libraries", library, "cachedStructures.pkl"), "rb") as f:
        cachedStructures_filtered = pickle.load(f)

    siriusDirectory = os.path.join("../data/libraries",library,"nf_output/fragmentationtrees/")

    if accepted_adduct is not None:
        matches_array = []
        for match in matches[1]:
            if count is not None:
                if len(matches_array) >= count:
                    break
            if data_dict_filtered[match[0]]['Adduct'] != accepted_adduct:
                continue
            if data_dict_filtered[match[1]]['Adduct'] != accepted_adduct:
                continue
            matches_array.append(match)
    else:
        matches_array = list(matches[1])
        if count is not None:
            matches_array = matches_array[:count]

    helpers = dict()
    for match in matches[1]:
        if match[0] not in helpers:
            helpers[match[0]] = []
        helpers[match[0]].append(match[1])

    return data_dict_filtered, cachedStructures_filtered, siriusDirectory, helpers, matches_array

def create_compound(main_id, args, data_dict_filtered, cachedStructures_filtered = {}, siriusDirectory = None, helpers = []):
    c = Compound(data_dict_filtered[main_id], cachedStructures_filtered.get(main_id, None), args)
    if siriusDirectory is not None:
        try:
            with open(os.path.join(siriusDirectory, main_id + "_fragmentationtree.json")) as f:
                sirius_json = json.load(f)
            c.apply_sirius(sirius_json)
        except:
            pass
    
    if len(helpers) > 0:
        for helper in helpers:
            helper_c = create_compound(helper, args, data_dict_filtered, cachedStructures_filtered, siriusDirectory, helpers=[])
            c.apply_helper(helper_c)
    return c

def create_modifiSiteLoc(main_id, modified_id, args, data_dict_filtered, cachedStructures_filtered = {}):
    main_compound = Compound(data_dict_filtered[main_id], cachedStructures_filtered.get(main_id, None), args)
    modified_compound = Compound(data_dict_filtered[modified_id], cachedStructures_filtered.get(modified_id, None), args)
    return ModificationSiteLocator(main_compound, modified_compound, args)

def multiprocessing_wrapper(func, library = "BERKELEY-LAB", accepted_adduct = None, count = None, processes = 16):
    data_dict_filtered, cachedStructures_filtered, siriusDirectory, helpers, matches_array = load(library, accepted_adduct, count)
    array = list(range(len(matches_array)))

    # add data_dict_filtered, cachedStructures_filtered, siriusDirectory, helpers, matches_array to shared memory
    par_func = partial(func, data_dict_filtered = data_dict_filtered, 
                       cachedStructures_filtered=cachedStructures_filtered,
                         siriusDirectory=siriusDirectory, helpers=helpers,
                           matches_array=matches_array)

    # Create a multiprocessing pool with desired number of processes
    pool = mp.Pool(processes=processes)  # Use 16 processes, adjust as needed

    results = []
    with tqdm(total=len(array), desc="Progress") as pbar:  # Initialize progress bar
        for result in pool.imap(par_func, array):
            results.append(result)
            pbar.update(1)
    results = filter(None, results)
    return results

def draw_result():
    pass