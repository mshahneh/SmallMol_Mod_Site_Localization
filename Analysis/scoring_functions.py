
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

args = {}

def compare_score_row(index, data_dict_filtered, cachedStructures_filtered, siriusDirectory, helpers, matches_array):
    """
    Compares the different scoring functions.
    """
    main_compound = create_compound(matches_array[index][1], args, data_dict_filtered, cachedStructures_filtered, siriusDirectory, helpers = [])
    modified_compound = create_compound(matches_array[index][0], args, data_dict_filtered, cachedStructures_filtered, siriusDirectory, helpers = [])
    site_locator = ModificationSiteLocator(main_compound, modified_compound, args)
    modification_site = utils.calculateModificationSites(modified_compound.structure, main_compound.structure, False)
    result = dict()
    for method in ["is_max", "dist_from_max", "average_dist_from_max", "average_dist"]:
        result["true_" + method] = site_locator.calculate_score(modification_site, method)
    
    # randomly set one of the atoms to be the modification site
    random_choice = np.zeros(len(main_compound.structure.GetAtoms()))
    random_choice[np.random.randint(0, len(random_choice))] = 1
    for method in ["is_max", "dist_from_max", "average_dist_from_max", "average_dist"]:
        result["random_choice" + method] = site_locator.calculate_score(modification_site, method, probabilities=random_choice)
    
    # random probabilities
    random_probabilities = np.random.rand(len(main_compound.structure.GetAtoms()))
    random_probabilities /= np.sum(random_probabilities)
    for method in ["is_max", "dist_from_max", "average_dist_from_max", "average_dist"]:
        result["random_probabilities" + method] = site_locator.calculate_score(modification_site, method, probabilities=random_probabilities)

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
    for method in ["is_max", "dist_from_max", "average_dist_from_max", "average_dist"]:
        # print the columns side by side
        print(df[["true_" + method, "random_choice" + method, "random_probabilities" + method]].describe())
    
    labels = ["is_max", "dist_from_max", "average_dist_from_max", "average_dist"]
    settings = ["true_", "random_choice", "random_probabilities"]
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

    
    
    