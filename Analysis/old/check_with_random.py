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



if __name__ == "__main__":
    libraries = {
    "GNPS-MSMLS": "https://external.gnps2.org/gnpslibrary/GNPS-MSMLS.json",
    # "GNPS-NIH-NATURALPRODUCTSLIBRARY_ROUND2_POSITIVE": "https://external.gnps2.org/gnpslibrary/GNPS-NIH-NATURALPRODUCTSLIBRARY_ROUND2_POSITIVE.json",
    # "GNPS-LIBRARY": "https://gnps-external.ucsd.edu/gnpslibrary/GNPS-LIBRARY.json"
    }

    for library in libraries:
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

        wb = xlsxwriter.Workbook("../data/libraries/" + library + "/results_random.xlsx")
        ws = wb.add_worksheet("results")
        df = pd.DataFrame(columns=["mol1ID", "mol2ID", "mol1smile", "mol2smile", "delta_mass", "#_matched_peaks", "#_shifted_peaks",
                                                  "#_unshifted_peaks", "score", "Closest_Max_Atom_Distance", "Count_Max", "Is_Max", "cosine", "random_guess", "random_prob", "url"])
        for match in tqdm(matches[1]):
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
                site = modSite.SiteLocator(molUsi, modifUsi, molSmiles)
                modifLoc = utils.calculateModificationSites(modifMol, molMol, False)
                peak_presence_only = True
                combine = True
                res = site.accuracy_score(modifLoc[0], peak_presence_only=peak_presence_only, combine=combine, return_all=True)
                # generate random probability array
                prob = np.zeros(site.molMol.GetNumAtoms())
                randInt = np.random.randint(0, site.molMol.GetNumAtoms())
                prob[randInt] = 1
                res2 = site.tempScore(modifLoc[0], prob, True)

                prb = np.random.rand(site.molMol.GetNumAtoms())
                prb = prb / prb.sum()
                res3 = site.tempScore(modifLoc[0], prb, True)
                df = pd.concat([df, pd.DataFrame.from_records([{"mol1ID": molUsi, "mol2ID": modifUsi, "mol1smile": molSmiles, "mol2smile": data_dict_filtered[m0]['Smiles'], 
                                                                "delta_mass": abs(float(data_dict_filtered[m0]['Precursor_MZ']) - float(data_dict_filtered[m1]['Precursor_MZ'])),
                                                                "#_matched_peaks": len(site.matchedPeaks), "#_shifted_peaks": len(site.shifted), "#_unshifted_peaks": len(site.unshifted),
                                                                "score": res['score'], "Closest_Max_Atom_Distance": res['closestMaxAtomDistance'],
                                                                "Count_Max": res['count'], "Is_Max": res['isMax'], "cosine":site.cosine, "random_guess":res2['score'],
                                                                "random_prob":res3['score'], "url":visualizer.make_url("http://reza.cs.ucr.edu/", molUsi, modifUsi, molSmiles, modifSmiles, args=None) }])], ignore_index=True)
            except:
                # print stack trace
                import traceback
                traceback.print_exc()
                pass

        for column in df.columns:
            ws.write(0, df.columns.get_loc(column), column)

        for index, row in df.iterrows():
            for column in df.columns:
                ws.write(index+1, df.columns.get_loc(column), row[column])
        wb.close()

        # write down the dataframe in picke format
        with open(library + "_run_scores.pkl", "wb") as f:
            pickle.dump(df, f)