
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


def compare_score_row(index, data_dict_filtered, cachedStructures_filtered, siriusDirectory, helpers, matches_array):
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
    peak_presence_only = True
    combine = True
    # calculate score 
    pre_helper = site.accuracy_score(modifLoc[0], peak_presence_only=peak_presence_only, combine=combine, return_all=True)
    try:
        molSirius = json.load(open(os.path.join(siriusDirectory, m1 + "_fragmentationtree.json")))
        site.apply_sirius(molSirius)
    except:
        print ("error finding sirius file for molecule")
        pass
    post_sirius = site.accuracy_score(modifLoc[0], peak_presence_only=peak_presence_only, combine=combine)
    scores_unshifted, scores_shifted = site.calculate_score(peak_presence_only)
    scores = site.distance_score(scores_unshifted, scores_shifted, combine)
    scores = np.array(scores)
    scores -= np.min(scores)
    if np.max(scores) != 0:
        scores /= np.max(scores)

    max_scores_unshifted, max_scores_shifted = site.calculate_score(peak_presence_only, consider_intensity=False, modificationSite=modifLoc[0])
    max_scores = site.distance_score(max_scores_unshifted, max_scores_shifted, combine)
    max_scores = np.array(max_scores)
    max_scores -= np.min(max_scores)
    if np.max(max_scores) != 0:
        max_scores /= np.max(max_scores)

    



    result = dict()
    for method in scoring_methods:
        result["true_" + method] = site_locator.calculate_score(modification_site, method, probabilities=scores)
    
    # randomly set one of the atoms to be the modification site
    random_choice = np.zeros(len(main_compound.structure.GetAtoms()))
    random_choice[np.random.randint(0, len(random_choice))] = 1
    for method in scoring_methods:
        result["random_choice" + method] = site_locator.calculate_score(modification_site, method, probabilities=random_choice)
    
    # random probabilities
    random_probabilities = np.random.rand(len(main_compound.structure.GetAtoms()))
    random_probabilities /= np.sum(random_probabilities)
    for method in scoring_methods:
        result["random_probabilities" + method] = site_locator.calculate_score(modification_site, method, probabilities=random_probabilities)
    
    for method in scoring_methods:
        result["best_" + method] = site_locator.calculate_score(modification_site, method, probabilities=max_scores)
    
    all_one = np.ones(len(main_compound.structure.GetAtoms()))
    all_one /= np.sum(all_one)
    for method in scoring_methods:
        result["all_one_" + method] = site_locator.calculate_score(modification_site, method, probabilities=all_one)

    return result 

if __name__ == '__main__':

    # if score_comparison.csv exists, load it
    try:
        df = pd.read_csv("score_comparison.csv")
    except:
    
        res = multiprocessing_wrapper(compare_score_row, library = "BERKELEY-LAB", accepted_adduct = "M+H", count = 50, processes = 8)
        
        # create pandas dataframe
        df = pd.DataFrame(res)
        df.to_csv("score_comparison.csv")

    print(df.head(5))

    # for each method, compare the statistics of the true and random scores
    for method in scoring_methods:
        # print the columns side by side
        print(df[["true_" + method, "random_choice" + method, "random_probabilities" + method, "best_" + method]].describe())
    
    labels = scoring_methods
    
    data = []
    for propability_setting in settings:
        row = []
        for method in labels:
            # calculate mean
            row.append(df[propability_setting + method].mean())
        
        data.append(row)
    
    # plot the data as spider radar plot and connect the points
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(111, polar=True)
    N = len(labels)
    theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
    for i in range(len(data)):
        ax.plot(theta, data[i], label=settings[i])
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_xticks(theta)
    ax.set_xticklabels(labels)
    ax.set_rlabel_position(0)
    ax.set_ylim(0, 1)
    ax.set_yticks(np.arange(0, 1, 0.1))
    ax.legend()
    plt.savefig("score_comparison.png")

    
    
    