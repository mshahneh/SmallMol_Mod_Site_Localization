
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

import visualizer as visualizer
import utils as utils
import fragmentation_py as fragmentation_py
import library_downloader as library_downloader
import SiteLocator as modSite
import handle_network as hn
import alignment as alignment

# library = "GNPS-LIBRARY"
# library = "GNPS-MSMLS"
# library = "GNPS-NIH-NATURALPRODUCTSLIBRARY_ROUND2_POSITIVE"
library = "BERKELEY-LAB" #: "https://external.gnps2.org/gnpslibrary/BERKELEY-LAB.json"

# if directory does not exist, create it
if not os.path.exists( os.path.join("../data/libraries",library)):
    url = "https://gnps-external.ucsd.edu/gnpslibrary/" + library + ".json"
    location = "../data/libraries/" + library + "/"
    library_downloader.download(url, location, 0.5, 0.1)

# load data_dict_filtered
with open(os.path.join("../data/libraries",library,"data_dict_filtered.pkl"), "rb") as f:
    data_dict_filtered = pickle.load(f)

# load matches
with open(os.path.join("../data/libraries",library,"matches.pkl"), "rb") as f:
    matches = pickle.load(f)

# load cachedStructures_filtered
with open(os.path.join("../data/libraries",library,"cachedStructures.pkl"), "rb") as f:
    cachedStructures_filtered = pickle.load(f)


def check_helper():
    helpers = dict()
    for match in matches[1]:
        if match[0] not in helpers:
            helpers[match[0]] = []
        helpers[match[0]].append(match[1])
    
    count = 0
    helperDirectory = os.path.join("../data/libraries",library,"nf_output/fragmentationtrees/")
    s_no_helper = 0
    s_helper = 0
    for i, match in tqdm(enumerate(list(matches[1]))):
        if count > 100:
            break
        try:
            m0, m1 = match
            if data_dict_filtered[m0]['Adduct'] != data_dict_filtered[m1]['Adduct'] or data_dict_filtered[m0]['Adduct'] != "M+H":
                continue
            molMol = cachedStructures_filtered[m1]
            modifMol = cachedStructures_filtered[m0]
            molUsi = hn.generate_usi(m1, library)
            modifUsi = hn.generate_usi(m0, library)
            molSmiles = data_dict_filtered[m1]['Smiles']
            modifSmiles = data_dict_filtered[m0]['Smiles']
            site = modSite.SiteLocator(data_dict_filtered[m1], data_dict_filtered[m0], molSmiles)
            modifLoc = utils.calculateModificationSites(modifMol, molMol, False)
            peak_presence_only = False
            combine = True
            consider_intensity = False
            # calculate score
            scores_unshifted, scores_shifted = modSite.calculate_score(peak_presence_only, consider_intensity)
            pre_helper = modSite.distance_score(scores_unshifted, scores_shifted, combine)

            pre_helper = site.accuracy_score(modifLoc[0], peak_presence_only=peak_presence_only, combine=combine, return_all=True)
            # print ("score is:", pre_helper['score'])
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
            post_helper = site.accuracy_score(modifLoc[0], peak_presence_only=peak_presence_only, combine=combine, return_all=True)
            if pre_helper['score'] != post_helper['score']:
                # print("Molecule: " + m1)
                # print("Modification: " + m0)
                # print("Pre-helper: " + str(pre_helper))
                # print("Post-helper: " + str(post_helper))
                # print("Helpers: " + str(helpers[m1]))
                # print("")
                s_no_helper += pre_helper['score']
                s_helper += post_helper['score']
                count += 1
        except:
            count += 1
            raise Exception("error")
            pass
    print(s_no_helper/count, s_helper/count)

check_helper()