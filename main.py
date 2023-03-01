import argparse
import SiteLocator as SiteLocator
import visualizer as visualizer
import utils as utils
import os
import json
import pickle
import requests
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Draw
from copy import deepcopy
import matplotlib.pyplot as plt
from IPython.display import SVG

def main(s1, s2, smiles, output, show_img=True):

    molMol = Chem.MolFromSmiles(smiles)

    modsite = SiteLocator.SiteLocator(s1, s2, molMol)
    scores_unshifted, scores_shifted = modsite.calculate_score()
    scores = modsite.distance_score(scores_unshifted, scores_shifted)
    if (max(scores.values()) > 0):
        scores = [scores[x] / max(scores.values()) for x in scores.keys()]
    else:
        scores = [0 for x in scores.keys()]
    
    if show_img:
        fig, ax = plt.subplots(1, figsize=(4,4))
        d2d = Draw.MolDraw2DSVG(250,200)
        colors = dict()
        for i in range(0, molMol.GetNumAtoms()):
            colors[i] = (1-(scores[i]**2), 1-(scores[i]**2), 1-(scores[i]**2))
        d2d.DrawMolecule(molMol, highlightAtoms=list(range(molMol.GetNumAtoms())), highlightAtomColors=colors)
        d2d.FinishDrawing()
        SVG(d2d.GetDrawingText())
    
    result = {}
    result["# matched peaks"] = len(modsite.matchedPeaks)
    result['# unshifted peaks'] = len(modsite.unshifted)
    result['# shifted peaks'] = len(modsite.shifted)
    result['image'] = d2d.GetDrawingText()

    visualizer.table_to_xlsx(result, output)


if __name__ == '__main__':
    # arguments
    parser = argparse.ArgumentParser(description='Substructure Assignment')
    parser.add_argument('--s1', type=str, default='mzspec:GNPS:TASK-e31d8bbfb26c4c3e8f159f06e37d65bb-spectra/specs_ms.mgf:scan:3189', help='gnps accepted usi for known molecule')
    parser.add_argument('--s2', type=str, default='mzspec:GNPS:TASK-e31d8bbfb26c4c3e8f159f06e37d65bb-spectra/specs_ms.mgf:scan:3399', help='gnps accepted usi for modified molecule')
    parser.add_argument('--smiles', type=str, default='OC(c1c(CCCCC)cc(OC(c2c(OC)cc(OC)cc2CCCCC)=O)cc1OC)=O', help='smiles for the known molecule')
    parser.add_argument('--output', type=str, default="output.xlsx", help='output file')
    parser.add_argument('--show_img', type=bool, default=True, help='show final image')
    args = parser.parse_args()

    main(args.s1, args.s2, args.smiles, args.output)