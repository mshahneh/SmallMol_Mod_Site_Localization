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
import copy

import os
import sys
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)

import visualizer as visualizer
import utils as utils
import handle_network as hn
import fragmentation_py as fragmentation_py
import library_downloader as library_downloader
import SiteLocator as modSite
import alignment as alignment
from Compound_n import Compound
from ModificationSiteLocator import ModificationSiteLocator

library ="BERKELEY-LAB"
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


def is_good(main, modified):
    molData = data_dict_filtered[main]
    modifiedData = data_dict_filtered[modified]
    molMol = cachedStructures_filtered[main]
    modifiedMol = cachedStructures_filtered[modified]

    if molData["Adduct"] != modifiedData["Adduct"] or molData["Adduct"] != "M+H":
        return -7

    if abs(float(molData["Precursor_MZ"]) - float(modifiedData["Precursor_MZ"])) < 20:
        return -1
    matchedPeaks = alignment.handle_alignment(molData, modifiedData)['matchedPeaks']
    if len(matchedPeaks) < 6:
        return -2
    
    if not modifiedMol.HasSubstructMatch(molMol):
        return -3
    
    main_c = Compound(molData, molMol)
    modified_c = Compound(modifiedData, modifiedMol)

    before = copy.deepcopy(main_c.peak_fragments_map)
    helperDirectory = os.path.join("../data/libraries", library, "nf_output/fragmentationtrees/")
    try:
        with open(os.path.join(helperDirectory, main + "_fragmentationtree.json")) as f:
            molSirius = json.load(f)
    except FileNotFoundError:
        return -4
    main_c.apply_sirius(molSirius)
    flag = False
    for i, peak in enumerate(main_c.peaks):
        if (len(main_c.peak_fragments_map[i]) == 1 and len(before[i]) != 1):
            flag = True
            break
    
    if flag:
        return -5
    
    site_locator = ModificationSiteLocator(main_c, modified_c)
    true_modif_loc = utils.calculateModificationSites(modifiedMol, molMol, False)
    score = site_locator.calculate_score(true_modif_loc, "temp")
    if score < 0.8:
        return -6
    
    return score

def get_good_matches(matches):
    good_matches = []
    diff = 0
    peak_count = 0
    not_sub = 0
    for match in tqdm(matches[1]):
        m0 = match[1]
        m1 = match[0]
        temp = is_good(m0, m1)
        if temp > 0:
            good_matches.append((temp, (m0, m1)))
        else:
            if temp == -1:
                diff += 1
            elif temp == -2:
                peak_count += 1
            elif temp == -3:
                not_sub += 1

    print (not_sub, peak_count, diff)
    return good_matches

good_matches = get_good_matches(matches)
good_matches.sort(reverse=True)
print(len(good_matches))
print(good_matches[0:10])
with open(os.path.join("./good_matches.pkl"), "wb") as f:
    pickle.dump(good_matches, f)