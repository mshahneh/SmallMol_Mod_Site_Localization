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

import visualizer as visualizer
import utils as utils
import handle_network as hn
import fragmentation_py as fragmentation_py
import library_downloader as library_downloader
import SiteLocator as modSite

import math
import multiprocessing as mp
from multiprocessing import Pool

libraries = {
    "GNPS-MSMLS": "https://external.gnps2.org/gnpslibrary/GNPS-MSMLS.json",
    "GNPS-NIH-NATURALPRODUCTSLIBRARY_ROUND2_POSITIVE": "https://external.gnps2.org/gnpslibrary/GNPS-NIH-NATURALPRODUCTSLIBRARY_ROUND2_POSITIVE.json",
    "GNPS-NIH-SMALLMOLECULEPHARMACOLOGICALLYACTIVE": "https://external.gnps2.org/gnpslibrary/GNPS-NIH-SMALLMOLECULEPHARMACOLOGICALLYACTIVE.json",
    "MIADB": "https://external.gnps2.org/gnpslibrary/MIADB.json",
    "BERKELEY-LAB": "https://external.gnps2.org/gnpslibrary/BERKELEY-LAB.json"
    # "GNPS-LIBRARY": "https://gnps-external.ucsd.edu/gnpslibrary/GNPS-LIBRARY.json"
}

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

helperDirectory = os.path.join("../data/libraries",library,"nf_output/fragmentationtrees/")
helpers = dict()
for match in matches[1]:
    if match[0] not in helpers:
        helpers[match[0]] = []
    helpers[match[0]].append(match[1])

print(len(helpers))
matches_array = list(matches[1])

# Define a function to create rows in the dataframe

# Define a function to be executed by each process
def process_element(element):
    try:
        m0, m1 = matches_array[element]
        if data_dict_filtered[m0]['Adduct'] != data_dict_filtered[m1]['Adduct'] or data_dict_filtered[m0]['Adduct'] != "M+H":
            return None
        molMol = cachedStructures_filtered[m1]
        modifMol = cachedStructures_filtered[m0]
        molUsi = hn.generate_usi(m1, library)
        modifUsi = hn.generate_usi(m0, library)
        molSmiles = data_dict_filtered[m1]['Smiles']
        modifSmiles = data_dict_filtered[m0]['Smiles']
        site = modSite.SiteLocator(data_dict_filtered[m1], data_dict_filtered[m0], molSmiles)
        modifLoc = utils.calculateModificationSites(modifMol, molMol, False)
        peak_presence_only = True
        combine = True
        # calculate score 
        pre_helper = site.accuracy_score(modifLoc[0], peak_presence_only=peak_presence_only, combine=combine, return_all=True)
        try:
            molSirius = json.load(open(os.path.join(helperDirectory, m1 + "_fragmentationtree.json")))
            site.apply_sirius(molSirius)
        except:
            print ("error finding sirius file for molecule")
            pass
        post_sirius = site.accuracy_score(modifLoc[0], peak_presence_only=peak_presence_only, combine=combine, return_all=True)
        for helper in helpers.get(m1, []):
            if helper != m0:
                helperFile = json.load(open(os.path.join(helperDirectory, helper + "_fragmentationtree.json")))
                try:
                    countUpdated = site.helper_molecule(data_dict_filtered[helper], data_dict_filtered[helper]['Smiles'], helperFile)
                    post_helper = site.accuracy_score(modifLoc[0], peak_presence_only=peak_presence_only, combine=combine, return_all=True)
                    # if pre_helper['score'] != post_helper['score']:
                    #     print(pre_helper['score'], post_helper['score'])
                except:
                    import traceback
                    traceback.print_exc()
                    pass
                break
        post_helper = site.accuracy_score(modifLoc[0], peak_presence_only=peak_presence_only, combine=combine, return_all=True)

        # generate random probability array 1-hot
        prob = np.zeros(site.molMol.GetNumAtoms())
        randInt = np.random.randint(0, site.molMol.GetNumAtoms())
        prob[randInt] = 1
        res2 = site.tempScore(modifLoc[0], prob, True)

        # generate random probability array distribution
        prb = np.random.rand(site.molMol.GetNumAtoms())
        prb = prb / prb.sum()
        res3 = site.tempScore(modifLoc[0], prb, True)

        # get max score
        maxScore = site.get_max_possible_score(modifLoc[0], peak_presence_only=peak_presence_only, combine=combine)
        
        row = {"mol1ID": molUsi, "mol2ID": modifUsi, "mol1smile": molSmiles, "mol2smile": data_dict_filtered[m0]['Smiles'], 
                                                        "delta_mass": abs(float(data_dict_filtered[m0]['Precursor_MZ']) - float(data_dict_filtered[m1]['Precursor_MZ'])),
                                                        "#_matched_peaks": len(site.matchedPeaks), "#_shifted_peaks": len(site.shifted), "#_unshifted_peaks": len(site.unshifted),
                                                        "Closest_Max_Atom_Distance": pre_helper['closestMaxAtomDistance'], "Count_Max": pre_helper['count'], "Is_Max": pre_helper['isMax'], "cosine":site.cosine, 
                                                        "pre_helper": float(pre_helper['score']), "post_sirius": float(post_sirius['score']) ,"post_helper": float(post_helper['score']), "best_score": maxScore, "random_guess":res2['score'], "random_prob":res3['score'], 
                                                        "url":visualizer.make_url("http://reza.cs.ucr.edu/", molUsi, modifUsi, molSmiles, modifSmiles, args=None) }
        return row
    except:
        # print stack trace
        import traceback
        traceback.print_exc()
        return None

if __name__ == '__main__':
    # Define your array of elements
    array = list(range(min(len(matches_array), 30000)))

    # Create a multiprocessing pool with desired number of processes
    pool = mp.Pool(processes=16)  # Use 16 processes, adjust as needed

    results = []
    with tqdm(total=len(array), desc="Progress") as pbar:  # Initialize progress bar
        for result in pool.imap(process_element, array):
            results.append(result)
            pbar.update(1)
    results = filter(None, results)

    # Close the pool to release resources
    pool.close()
    pool.join()

    # Create a dataframe from the results
    df = pd.DataFrame(results)
    resultColumns = ['pre_helper', 'post_sirius','post_helper', "best_score", "random_guess", "random_prob"]
    print(df[resultColumns].describe())

    # select the rows that have a score difference
    df2 = df[df['pre_helper'] != df['post_helper']]
    print(df2[resultColumns].describe())
    
    df_stats = df[resultColumns].describe()
    bar_plot = df_stats.loc[['std', 'mean', '25%', '50%', '75%']].plot(kind='bar', legend=True)
    plt.title('Bar Plot of Descriptive Statistics')
    plt.xlabel('Statistics')
    plt.ylabel('Value')
    plt.legend(loc='upper right')
    plt.savefig('df_stats_box_plot.png', bbox_inches='tight')
    plt.close()

    df2_stats = df2[resultColumns].describe()
    bar_plot = df2_stats.loc[['std', 'mean', '25%', '50%', '75%']].plot(kind='bar', legend=True)
    plt.title('Bar Plot of Descriptive Statistics')
    plt.xlabel('Statistics')
    plt.ylabel('Value')
    plt.legend(loc='upper right')
    plt.savefig('df2_stats_box_plot.png', bbox_inches='tight')
    plt.close()
