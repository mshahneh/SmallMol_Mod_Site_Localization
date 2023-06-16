
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

args = {}
scoring_methods = ["is_max", "dist_from_max", "average_dist_from_max", "average_dist", "temp"]
settings = ["true_", "random_choice", "random_probabilities", "best_", "all_one_"]


def peak_ambiguity_row(index, data_dict_filtered, cachedStructures_filtered, siriusDirectory, helpers, matches_array):
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
    pre_helper = site.accuracy_score(modifLoc[0], peak_presence_only=peak_presence_only, combine=combine, return_all=True)
    try:
        molSirius = json.load(open(os.path.join(siriusDirectory, m1 + "_fragmentationtree.json")))
        site.apply_sirius(molSirius)
    except:
        print ("error finding sirius file for molecule")
        pass
    post_sirius = site.accuracy_score(modifLoc[0], peak_presence_only=peak_presence_only, combine=combine)

    average_ambiguity = 0
    for peak in range(len(main_compound.peak_fragments_map)):
        average_ambiguity += len(main_compound.peak_fragments_map[peak])
    if len(main_compound.peak_fragments_map) != 0:  
        average_ambiguity /= len(main_compound.peak_fragments_map)
    else:
        return None

    result = {"ambiguity": average_ambiguity, "score": post_sirius}

    return result 

if __name__ == '__main__':

    # if peak_ambiguity_effect.csv exists, load it
    try:
        df = pd.read_csv("peak_ambiguity_effect.csv")
    except:
    
        res = multiprocessing_wrapper(peak_ambiguity_row, library = "BERKELEY-LAB", accepted_adduct = "M+H", count = 500, processes = 8)
        
        # create pandas dataframe
        df = pd.DataFrame(res)
        df.to_csv("peak_ambiguity_effect.csv")

    print(df.head(5))

    # drop score = 0
    df = df[df["score"] != 0]

    # group by ambiguity with 0.5 range and calculate mean
    # df = df.groupby(pd.cut(df["ambiguity"], np.arange(0, 10, 0.5))).mean()

    scores = df["score"].tolist()
    ambiguities = df["ambiguity"].tolist()

    
    plt.figure(figsize=(40, 40))

    # plot reverse log x axis
    plt.xscale('function', functions=(lambda x: np.log10(x+0.001), lambda x: 10**x))

    # plot data
    plt.scatter(ambiguities, scores, s=190)

    # change size of ticks and labels
    plt.tick_params(axis='both', which='major', labelsize=100)
    plt.tick_params(axis='both', which='minor', labelsize=100)

    # set ticks
    plt.xticks([1, 2, 4, 8, 16, 32, 64], ["1", "2", "4", "8", "16", "32", "64"])

    # change size of axis labels
    plt.xlabel("Average peak ambiguity", fontsize=140)
    plt.ylabel("Weighted Average Score", fontsize=140)


    # plt.figure(figsize=(10, 10))
    # plt.plot(df["ambiguity"], df["score"])
    # plt.xlabel("Average peak ambiguity")
    # plt.ylabel("Score")

   
    plt.savefig("peak_ambiguity_effect.png")