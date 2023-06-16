
from base_tester import multiprocessing_wrapper, create_compound
from ModificationSiteLocator import ModificationSiteLocator
import utils as utils
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, RegularPolygon
from matplotlib.path import Path
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D
import SiteLocator as modSite
import json
import os

args={}

def sirius_effect_row(index, data_dict_filtered, cachedStructures_filtered, siriusDirectory, helpers, matches_array):
    """
    Compares the different scoring functions.
    """
    main_compound = create_compound(matches_array[index][1], args, data_dict_filtered, cachedStructures_filtered, siriusDirectory, helpers = [])
    modified_compound = create_compound(matches_array[index][0], args, data_dict_filtered, cachedStructures_filtered, siriusDirectory, helpers = [])
    site_locator = ModificationSiteLocator(main_compound, modified_compound, args)
    modification_site = utils.calculateModificationSites(modified_compound.structure, main_compound.structure, False)[0]
    
    m0, m1 = matches_array[index]
    molMol = cachedStructures_filtered[m1]
    modifMol = cachedStructures_filtered[m0]
    molSmiles = data_dict_filtered[m1]['Smiles']
    modifSmiles = data_dict_filtered[m0]['Smiles']
    site = modSite.SiteLocator(data_dict_filtered[m1], data_dict_filtered[m0], molSmiles)
    modifLoc = utils.calculateModificationSites(modifMol, molMol, False)
    peak_presence_only = False
    combine = False
    # calculate score 
    post_sitius_probs = site_locator.generate_probabilities(shifted_only=(not combine), PPO=peak_presence_only, CI=False, modificationSite=modifLoc[0])
    post_sirius = site_locator.calculate_score(modification_site, method, probabilities=post_sitius_probs)
    pre_helper = site.accuracy_score(modifLoc[0], peak_presence_only=peak_presence_only, combine=combine, return_all=True)
    try:
        molSirius = json.load(open(os.path.join(siriusDirectory, m1 + "_fragmentationtree.json")))
        site.apply_sirius(molSirius)
    except:
        print ("error finding sirius file for molecule")
        pass
    post_sitius_probs = site_locator.generate_probabilities(shifted_only=(not combine), PPO=peak_presence_only, CI=False, modificationSite=modifLoc[0])
    post_sirius = site_locator.calculate_score(modification_site, method, probabilities=post_sitius_probs)

    result = {"pre": pre_helper, "post": post_sirius}

    random_choice = np.zeros(len(main_compound.structure.GetAtoms()))
    random_choice[np.random.randint(0, len(random_choice))] = 1
    result["random_choice" + method] = site_locator.calculate_score(modification_site, method, probabilities=random_choice)
    
    # random probabilities
    random_probabilities = np.random.rand(len(main_compound.structure.GetAtoms()))
    random_probabilities /= np.sum(random_probabilities)
    # for method in scoring_methods:
    #     result["random_probabilities" + method] = site_locator.calculate_score(modification_site, method, probabilities=random_probabilities)

    return result 

if __name__ == '__main__':

    # if peak_ambiguity_effect.csv exists, load it
    try:
        df = pd.read_csv("sirius_effect.csv")
    except:
    
        res = multiprocessing_wrapper(sirius_effect_row, library = "BERKELEY-LAB", accepted_adduct = "M+H", count = 2500, processes = 8)
        
        # create pandas dataframe
        df = pd.DataFrame(res)
        df.to_csv("sirius_effect.csv")
