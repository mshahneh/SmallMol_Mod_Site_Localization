from dash import Dash, html, dcc, Input, Output, State, dash_table
from dash.exceptions import PreventUpdate
import plotly.graph_objects as go
import plotly.express as px
from dash.dash_table.Format import Format, Scheme, Trim
from urllib.request import urlopen
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from tqdm import tqdm
import pandas as pd
import subprocess
import traceback
import pickle
import json
import os
import base64
import visualizer as visualizer
import utils as utils
import fragmentation_py as fragmentation_py
import SiteLocator as modSite
from IPython.display import SVG


df = pd.DataFrame(columns=["mol1ID", "mol1Weight", "mol1Library", "mol2ID", "mol2Weight", "mol1Library", "score", "# matched peaks", "# shifted peaks", "# unshifted peaks"])
columns=[{"name": i, "id": i} for i in df.columns]
columns.append({"name": "structure", "id": "structure", "presentation": "markdown"})
columns[2]["format"] = Format(precision=2, scheme=Scheme.fixed)
# print(columns)

library = "GNPS-MSMLS"
with open('data/libraries/{}/data_dict_filtered.pkl'.format(library), 'rb') as f:
    data_dict_filtered = pickle.load(f)

with open('data/libraries/{}/matches.pkl'.format(library), 'rb') as f:
    matches = pickle.load(f)

with open('data/libraries/{}/cachedStructures.pkl'.format(library), 'rb') as f:
    cachedStructures = pickle.load(f)
app = Dash()

siteLocator = None

app.layout = html.Div(id = 'parent', children = [
        html.Div(id = 'inputs', children = [
                    dcc.Input(
            id="USI1",
            type="text",
            placeholder="USI1",
            value = "mzspec:GNPS:TASK-5700dee92610412ea452a4262add2b93-f.MSV000086107/ccms_peak/VVP3-2_EtAc_MeOh.mzML:scan:2303",
        ),
                    dcc.Input(
            id="USI2",
            type="text",
            placeholder="USI2 (modified bigger molecule)",
            value = "mzspec:GNPS:TASK-5700dee92610412ea452a4262add2b93-f.MSV000086107/ccms_peak/VVP77_MeOh.mzML:scan:3024",
        ),
                    dcc.Input(
            id="SMILES1",
            type="text",
            placeholder="SMILES1",
            value = "OC1=CC=CC=C1C2=NC(C3N(C)C(C(O)=O)CS3)CS2",
        ),
                    dcc.Input(
            id="SMILES2",
            type="text",
            placeholder="SMILES2 (optional)",
            value = "OC1=CC=CC=C1C2=NC(C3N(C)C(C(OC)=O)CS3)CS2",
        ),
        ]),
        html.Div(id = 'retrive', children = [
            html.Button('Refresh', id='refresh', n_clicks=0),
            html.P(id='retrive_status'),
        ]),
        html.Img(id = 'prediction'),
        html.Img(id = 'substr'),
        #bar chart
        dcc.Graph(id='peaks'),
        # text
        html.Div(id='debug_text'),
        dcc.Store(id='siteLocatorObj'),
    ])

@app.callback(
    [
        Output('substr', 'src'),
        Output('prediction', 'src'),
        Output('peaks', 'figure'),
        Output('retrive_status', 'children'),
        Output('siteLocatorObj', 'data'),
    ],
    Input('refresh', 'n_clicks'),
    State('USI1', 'value'),
    State('USI2', 'value'),
    State('SMILES1', 'value'),
    State('SMILES2', 'value'),
    )
def update_output(n_clicks, usi1, usi2, smiles1, smiles2):
    if (n_clicks == 0 or n_clicks is None):
        raise PreventUpdate
    
    if (usi1 is None or usi2 is None or smiles1 is None):
        return None, None, None, "Please enter USI1, USI2, and SMILES1.", None
    
    mol1 = Chem.MolFromSmiles(smiles1)
    siteLocator = modSite.SiteLocator(usi1, usi2, smiles1)
    if siteLocator.molMeta['precursor_mz'] > siteLocator.modifMeta['precursor_mz']:
        return None, None, None, "Not supported, Modified molecule is smaller than the original molecule."
    scores_unshifted, scores_shifted = siteLocator.calculate_score()
    scores = siteLocator.distance_score(scores_unshifted, scores_shifted)
    if (max(scores.values()) > 0):
        scores = [scores[x] / max(scores.values()) for x in scores.keys()]
    else:
        scores = [0 for x in scores.keys()]
    
    if (smiles2 is not None and len(smiles2) > 0):
        svg1 = visualizer.molToSVG(Chem.MolFromSmiles(smiles2), mol1, True)
    else:
        svg1 = None
    
    svg2 = visualizer.highlightScores(mol1, scores)

    fig = go.Figure()
    x1 = []
    y1 = []
    for peak in siteLocator.molPeaks:
        x1.append(peak[0])
        y1.append(peak[1])
    

    topPeakCount = 50
    # get index of top peaks in y1
    topPeaksInx = sorted(range(len(y1)), key=lambda i: y1[i])[-topPeakCount:]
    types = ['matched_shifted', 'matched_unshifted', 'unmatched']
    typesInx = {'matched_shifted': [], 'matched_unshifted': [], 'unmatched': []}
    for i in topPeaksInx:
        flag = False
        for j in siteLocator.matchedPeaks:
            if (j[0] == i):
                if (abs(siteLocator.molPeaks[i][0] - siteLocator.modifPeaks[j[1]][0]) > 0.1):
                    typesInx['matched_shifted'].append(i)
                else:
                    typesInx['matched_unshifted'].append(i)
                flag = True
                break
        if (not flag):
            typesInx['unmatched'].append(i)
    
    for i in typesInx:
        x1_ = [round(x1[j], 2) for j in typesInx[i]]
        y1_ = [y1[j] for j in typesInx[i]]
        y1_ = [x / max(y1_)*100 for x in y1_]
        fig.add_trace(go.Bar(x=x1_, y=y1_, name=i, width=2))

    x2 = []
    y2 = []
    for peak in siteLocator.modifPeaks:
        x2.append(peak[0])
        y2.append(peak[1])

    # get index of top 10 peaks in y2
    topPeaksInx = sorted(range(len(y2)), key=lambda i: y2[i])[-topPeakCount:]
    x2 = [round(x2[i], 2) for i in topPeaksInx]
    y2 = [y2[i] for i in topPeaksInx]
    y2 = [-x / max(y2)*100 for x in y2]

    fig.add_trace(go.Bar(x=x2, y=y2, name='modified molecule peaks', width=2))
    fig.update_layout(
        title="Peaks",
        bargap=0,
        xaxis_title="m/z",
        yaxis_title="intensity",
        legend_title="Peak Type",
    )

    return dash_svg(svg1), dash_svg(svg2), fig, "Success", base64.b64encode(pickle.dumps(siteLocator)).decode()

## update debugtext if click on bar chart
@app.callback(
    Output('debug_text', 'children'),
    Input('peaks', 'clickData'),
    State('siteLocatorObj', 'data'), prevent_initial_call=True)
def display_click_data(clickData, siteLocatorObj):
    if clickData:
        try:
            siteLocator = pickle.loads(base64.b64decode(siteLocatorObj))
            structures = siteLocator.get_structures_per_peak(float(clickData['points'][0]['x']))
            res = []
            for structure in structures:
                svg = dash_svg(visualizer.molToSVG(siteLocator.molMol, Chem.MolFromSmiles(structure)))
                res.append(html.Img(src=svg, style={'width': '300px'}))
            return res
                
        except:
            return "siteLocator object not found"
    return None
   

def dash_svg(text):
    """
    Generates a svg link for dash from a avg text.
    """
    try:
        svg = base64.b64encode(text.encode('utf-8'))
        return 'data:image/svg+xml;base64,{}'.format(svg.decode())
    except:
        return None

if __name__ == '__main__':

    app.run_server(debug=True )#port=80, host="0.0.0.0")