  
from prettytable import PrettyTable
from urllib.request import urlopen
import matplotlib.pyplot as plt
from rdkit import Chem
from tqdm import tqdm
import pandas as pd
import pickle
import utils
import json
import os

def disable_rdkit_logging():
    """
    Disables RDKit whiny logging.
    """
    import rdkit.rdBase as rkrb
    import rdkit.RDLogger as rkl
    logger = rkl.logger()
    logger.setLevel(rkl.ERROR)
    rkrb.DisableLog('rdApp.error')

def get_gnps_library(url):
    """
    Downloads the GNPS library from the GNPS website.
    """
    print("Downloading GNPS library from: " + url)
    response = urlopen(url)
    data_json = json.loads(response.read())
    print ("Number of compounds in the library: " + str(len(data_json)))
    return data_json

def json_to_dict(data_json):
    """
    Converts the GNPS library JSON to a dictionary.
    """
    data_dict = dict()
    for i in range(len(data_json)):
        data_dict[data_json[i]['spectrum_id']] = data_json[i]
    return data_dict

def calculate_matches(data_dict, weight_threshold = 500, difference_threshold_rate = 0.5):
    """
    Calculates the matches between the library.
    It is assumed that the library is passed as a dictionary with unique key identifiers.
    Input:
        data_dict: dictionary with the library
        weight_threshold: maximum weight to consider
        difference_threshold_rate: maximum difference between two weights to consider
    Output:
    """
    matches = {}
    cachedStructures = dict()
    data_ids = list(data_dict.keys())
    weight_threshold = 500
    difference_threshold_rate = 0.5
    for i in tqdm(range(len(data_ids))):   
        compound1 = data_ids[i]
        # print(compound1)
        # if compound1 == "CCMSLIB00005884654" or compound1 == "CCMSLIB00005883712":
        #     print("here!", i)
        try:
            w1 = float(data_dict[compound1]['Precursor_MZ'])
        except:
            continue
        # if compound1 == "CCMSLIB00005884654" or compound1 == "CCMSLIB00005883712":
        #     print("here2!", w1, weight_threshold, data_dict[compound1]['Smiles'], )
        if w1 > weight_threshold or data_dict[compound1]['Smiles'] == 'N/A' or data_dict[compound1]['Smiles'] == ' ':
            continue
        
        for j in range(len(data_ids)):

            # if compound1 == "CCMSLIB00005884654" or compound1 == "CCMSLIB00005883712":
            #     print(".", end="")
            compound2 = data_ids[j]
            # if j == 9840:
            #     print("WHAT!?")
            # if compound1 == "CCMSLIB00005884654" and compound2 == "CCMSLIB00005883712":
            #     print("here!")
            # if compound1 == "CCMSLIB00005883712" and compound2 == "CCMSLIB00005884654":
            #     print("here!")
            try:
                w2 = float(data_dict[compound2]['Precursor_MZ'])
            except:
                continue
            if w2 > weight_threshold or data_dict[compound2]['Smiles'] == 'N/A' or data_dict[compound2]['Smiles'] == ' ':
                continue
            # if compound1 == "CCMSLIB00005883712" and compound2 == "CCMSLIB00005884654":
            #     print("here!")
            #     print(w1, w2, difference_threshold_rate*w2)
            
            if abs(w1 - w2) > difference_threshold_rate * min(w1, w2):
                continue

            try:
                if compound1 not in cachedStructures:
                    cachedStructures[compound1] = Chem.MolFromSmiles(data_dict[compound1]['Smiles'], sanitize=False)
                if compound2 not in cachedStructures:
                    cachedStructures[compound2] = Chem.MolFromSmiles(data_dict[compound2]['Smiles'], sanitize=False)

                m1 = cachedStructures[compound1]
                m2 = cachedStructures[compound2]
                
                # if compound1 == "CCMSLIB00005883712" and compound2 == "CCMSLIB00005884654":
                #     print("here")
                #     print (m1.GetNumAtoms())
                #     print (m2.GetNumAtoms())
                #     print (m1.HasSubstructMatch(m2))

                if m1.GetNumAtoms() > m2.GetNumAtoms() and m1.HasSubstructMatch(m2):
                    numModificationSites = len(utils.calculateModificationSites(m1, m2))
                    if numModificationSites not in matches:
                        matches[numModificationSites] = set()
                    matches[numModificationSites].add((compound1, compound2))

            except:
                pass
    
    return matches, cachedStructures

def main(url, output, weight_threshold, difference_threshold_rate):
    """
    Main function.
    """
    disable_rdkit_logging()
    data_json = get_gnps_library(url)
    data_dict = json_to_dict(data_json)
    matches, cachedStructures = calculate_matches(data_dict, weight_threshold, difference_threshold_rate)
    
    # get compound ids that are in matches
    unique_ids = set()
    for key in matches:
        for match in matches[key]:
            unique_ids.add(match[0])
            unique_ids.add(match[1])
    
    # filter data_dict to only include unique_ids
    data_dict_filtered = dict()
    for key in data_dict:
        if key in unique_ids:
            data_dict_filtered[key] = data_dict[key]

    print("Number of compounds in the library: " + str(len(data_dict_filtered)))
    t = PrettyTable(['# modification sites', '# matches'])
    for key in matches:
        t.add_row([key, len(matches[key])])
    print(t)
    print("Number of unique ids in matches: " + str(len(unique_ids)))
    
    # save data_dict_filtered
    with open(os.path.join(output, "data_dict_filtered.pkl"), "wb") as f:
        pickle.dump(data_dict_filtered, f)
    
    # save matches
    with open(os.path.join(output, "matches.pkl"), "wb") as f:
        pickle.dump(matches, f)

    # save cachedStructures
    with open(os.path.join(output, "cachedStructures.pkl"), "wb") as f:
        pickle.dump(cachedStructures, f)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='GNPS library downloader')
    parser.add_argument('--url', type=str, default="https://gnps-external.ucsd.edu/gnpslibrary/GNPS-LIBRARY.json", help='URL to download the JSON library from')
    parser.add_argument('--output', type=str, default="Data", help='Output path')
    parser.add_argument('--weight_threshold', type=float, default=500, help='Maximum weight to consider')
    parser.add_argument('--difference_threshold_rate', type=float, default=0.5, help='Maximum difference between two weights to consider')
    args = parser.parse_args()
    main(**vars(args))