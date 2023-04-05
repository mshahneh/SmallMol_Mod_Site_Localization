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
from furl import furl
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
        dcc.Location(id='url', refresh=False),
        html.Div(id = 'inputs', children = [
            html.Div(children = [
                dcc.Input(
                id="USI1",
                type="text",
                placeholder="USI1",
                value = "mzspec:GNPS:TASK-5700dee92610412ea452a4262add2b93-f.MSV000086107/ccms_peak/VVP3-2_EtAc_MeOh.mzML:scan:2303",
                style={'width':'23%'}
            ),
                        dcc.Input(
                id="USI2",
                type="text",
                placeholder="USI2 (modified bigger molecule)",
                value = "mzspec:GNPS:TASK-5700dee92610412ea452a4262add2b93-f.MSV000086107/ccms_peak/VVP77_MeOh.mzML:scan:3024",
                style={'width':'23%'}
            ),
                        dcc.Input(
                id="SMILES1",
                type="text",
                placeholder="SMILES1",
                value = "OC1=CC=CC=C1C2=NC(C3N(C)C(C(O)=O)CS3)CS2",
                style={'width':'23%'}
            ),
                        dcc.Input(
                id="SMILES2",
                type="text",
                placeholder="SMILES2 (optional)",
                value = "OC1=CC=CC=C1C2=NC(C3N(C)C(C(OC)=O)CS3)CS2",
                style={'width':'23%'}
            ),], style = {'display': 'flex', 'flex-direction': 'row', 'flex-wrap': 'wrap', 'justify-content': 'space-around', 'width': '100%', 'min-height': '5vh'}),
        
            html.Div(children = [
                dcc.Checklist(
                        id='options',
                        options=[
                            {'label': 'presense only', 'value': 'presense'},
                            {'label': 'subtract unshifted peaks', 'value': 'combine'},
                            {'label': 'Intensity in score', 'value': 'intensity'},
                    ], labelStyle={'display': 'inline-block', 'margin-right': '1vw'}),
                html.Div(id = 'arguments', children = [
                    dcc.Dropdown(['intensity', 'top_k', 'none'], 'top_k', id='filter_peaks_method'),
                    dcc.Input(
                        id="filter_peaks_variable",
                        type="text",
                        placeholder="filter peaks variable",
                        value = "50",
                        style={'margin-right': '1vw'}
                    ),
                    dcc.Input(
                        id="mz_tolerance",
                        type="text",
                        placeholder="mz_tolerance",
                        value = "0.1",
                    ),
                        dcc.Input(
                        id="ppm",
                        type="text",
                        placeholder="ppm",
                        value = "1.01",
                    )
                ], style = {'display': 'flex', 'flex-direction': 'row', 'justify-content': 'center', 'min-height': '5vh', 'align-items': 'center', 'margin-top': '1vh'}),
            ], style = {'display': 'flex', 'flex-direction': 'row', 'justify-content': 'center', 'align-items': 'center'}),
        ]),

        html.Div(id = 'retrive', children = [
            html.Button('Refresh', id='refresh', n_clicks=0),
            html.P(id='retrive_status'),
        ], style = {'display': 'flex', 'flex-direction': 'row', 'justify-content': 'center', 'align-items': 'center'}),
        
        html.Div(children = [
            html.Div(id = 'stats', style={"margin-right": "1vw"}),
            html.Div(children = [
                html.Img(id = 'prediction'),
                html.Img(id = 'substr'),
            ], style = {'display': 'flex', 'flex-direction': 'row', 'justify-content': 'center', 'align-items': 'center'})
        ],  style = {'display': 'flex', 'flex-direction': 'row', 'justify-content': 'center', 'align-items': 'center'}),

        html.Div(id = "peaks-component", children = [
            dcc.Graph(id='peaks', style = {'width': '100%'}),
            html.Div(children = [
                dcc.Slider(
                    id='peak_filter_slider',
                    min=0,
                    max=100,
                    step=1,
                    value=50,
                    marks=None,
                    tooltip={"placement": "bottom", "always_visible": True},
                )
            ], style = {'width': '30%', 'height': '5vh'}),
        ], style = {'display': 'none'}),
        # text
        html.Div(id='peak_info'),
        dcc.Store(id='siteLocatorObj'),
    ])

@app.callback(
    [Output('peaks', 'figure'),
     Output('peaks-component', 'style')],
    Input('siteLocatorObj', 'data'),
    Input('peak_filter_slider', 'value'), prevent_initial_call=True)
def update_peaks(data, value):
    print("update_peaks")
    if (data == None):
        return {}, {'display': 'none'}
    siteLocator = pickle.loads(base64.b64decode(data))
    fig = go.Figure()
    typesInx = {'matched_shifted': [], 'matched_unshifted': [], 'unmatched': []}
    x1 = []
    y1 = []
    for peak in siteLocator.molPeaks:
        x1.append(peak[0])
        y1.append(peak[1])
    

    topPeakCount = value
    topPeaksInxModif = sorted(range(len(y1)), key=lambda i: y1[i])[-topPeakCount:]
    for i in topPeaksInxModif:
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
        y1_ = [x / max(y1_) * 100 for x in y1_]
        indicis = typesInx[i]
        if (i == 'unmatched'):
            fig.add_trace(go.Bar(x=x1_, y=y1_, hovertext=indicis, name=i, width=1, visible='legendonly'))
        else:
            fig.add_trace(go.Bar(x=x1_, y=y1_, hovertext=indicis, name=i, width=1))

    typesInx = {'matched_shifted': [], 'matched_unshifted': [], 'unmatched': []}
    x2 = []
    y2 = []
    for peak in siteLocator.modifPeaks:
        x2.append(peak[0])
        y2.append(peak[1])
    
    topPeaksInxModif = sorted(range(len(y2)), key=lambda i: y2[i])[-topPeakCount:]
    for i in topPeaksInxModif:
        flag = False
        for j in siteLocator.matchedPeaks:
            if (j[1] == i):
                if (abs(siteLocator.molPeaks[j[0]][0] - siteLocator.modifPeaks[j[1]][0]) > 0.1):
                    typesInx['matched_shifted'].append(i)
                else:
                    typesInx['matched_unshifted'].append(i)
                flag = True
                break
        if (not flag):
            typesInx['unmatched'].append(i)
    

    for i in typesInx:
        x1_ = [round(x2[j], 2) for j in typesInx[i]]
        y2_ = [y2[j] for j in typesInx[i]]
        y2_ = [-j / max(y2_) * 100 for j in y2_]
        indicis = typesInx[i]
        if (i == 'unmatched'):
            fig.add_trace(go.Bar(x=x1_, y=y2_, hovertext=indicis, name=i, width=1, visible='legendonly'))
        else:
            fig.add_trace(go.Bar(x=x1_, y=y2_, hovertext=indicis, name=i, width=1))
    
    fig.update_layout(
        title="Peaks",
        bargap=0,
        xaxis_title="m/z",
        yaxis_title="intensity",
        legend_title="Peak Type",
    )
    return fig, {'display': 'flex', 'flex-direction': 'column', 'justify-content': 'center', 'width': '100%', 'align-items': 'center', 'margin-top': '-4vh'}


@app.callback(
    Output('url', 'search'),
    Input('refresh', 'n_clicks'),
    State('USI1', 'value'),
    State('USI2', 'value'),
    State('SMILES1', 'value'),
    State('SMILES2', 'value'),
    State('arguments', 'children'),
    State('options', 'value'),
    prevent_initial_call=True)
def update_url(n_clicks, usi1, usi2, smiles1, smiles2, arguments, options):
    if (n_clicks == 0 or n_clicks is None):
        raise PreventUpdate

    print("update_url")
    filter_peaks_method = arguments[0]['props']['value']
    filter_peaks_variable = float(arguments[1]['props']['value'])
    mz_tolerance = float(arguments[2]['props']['value'])
    ppm = float(arguments[3]['props']['value'])
    
    # generate url for the current state
    url = '?USI1=' + usi1 + '&USI2=' + usi2 + '&SMILES1=' + smiles1 + '&SMILES2=' + smiles2 + '&filter_peaks_method=' + filter_peaks_method + '&filter_peaks_variable=' + str(filter_peaks_variable) + '&mz_tolerance=' + str(mz_tolerance) + '&ppm=' + str(ppm) + '&options=' + str(options)
    return url


@app.callback(
    [
        Output('substr', 'src'),
        Output('prediction', 'src'),
        Output('retrive_status', 'children'),
        Output('siteLocatorObj', 'data'),
        Output('stats', 'children'),
    ],
    Input('refresh', 'n_clicks'),
    State('USI1', 'value'),
    State('USI2', 'value'),
    State('SMILES1', 'value'),
    State('SMILES2', 'value'),
    State('arguments', 'children'),
    State('options', 'value')
    )
def update_output(n_clicks, usi1, usi2, smiles1, smiles2, arguments, options):
    if (n_clicks == 0 or n_clicks is None):
        raise PreventUpdate
    print("update_output")
    filter_peaks_method = arguments[0]['props']['value']
    filter_peaks_variable = float(arguments[1]['props']['value'])
    mz_tolerance = float(arguments[2]['props']['value'])
    ppm = float(arguments[3]['props']['value'])
    presense = False
    combine = False
    consider_intensity = False
    if options is not None:
        if 'presense' in options:
            presense = True
        if 'combine' in options:
            combine = True
        if 'intensity' in options:
            consider_intensity = True

    print(filter_peaks_method, filter_peaks_variable, mz_tolerance, ppm)

    if (n_clicks == 0 or n_clicks is None):
        raise PreventUpdate
    
    if (usi1 is None or usi2 is None or smiles1 is None):
        return None, None, "Please enter USI1, USI2, and SMILES1.", None, None
    
    mol1 = Chem.MolFromSmiles(smiles1)

    args = {
        'filter_peaks_method': filter_peaks_method,
        'filter_peaks_variable': filter_peaks_variable,
        'mz_tolerance': mz_tolerance,
        'ppm': ppm,
    }

    siteLocator = modSite.SiteLocator(usi1, usi2, mol1, args)
    if siteLocator.molPrecursorMz > siteLocator.modifPrecursorMz:
        return None, None, None, "Not supported, Modified molecule is smaller than the original molecule."
    scores_unshifted, scores_shifted = siteLocator.calculate_score(peak_presence_only = presense, consider_intensity = consider_intensity)
    scores = siteLocator.distance_score(scores_unshifted, scores_shifted, combine = combine)
    
    isMax = None
    score = None
    if (smiles2 is not None and len(smiles2) > 0):
        svg1 = visualizer.molToSVG(Chem.MolFromSmiles(smiles2), mol1, True)
        modifLoc = list(utils.calculateModificationSites(Chem.MolFromSmiles(smiles2), mol1, False))
        accuracy_score = siteLocator.accuracy_score(modifLoc[0], peak_presence_only=presense, combine=combine, return_all=True, consider_intensity=consider_intensity)
        isMax = accuracy_score['isMax']
        score = accuracy_score['score']
        stats =  html.Div([
            html.P("cosine: " + str(round(siteLocator.cosine, 3)), style = {'margin-left': '2vw'}),
            html.P("Score: " + str(round(score, 3)), style = {'margin-left': '2vw'}),
            html.P("is Max: " + str(isMax), style = {'margin-left': '2vw'}),
            html.P("#matches: " + str(len(siteLocator.matchedPeaks)), style = {'margin-left': '2vw'}),
            html.P("#shifted: " + str(len(siteLocator.shifted)), style = {'margin-left': '2vw'}),
        ], style = {'display': 'flex', 'flex-direction': 'row', 'flex-wrap': 'wrap', 'justify-content': 'left', 'width': '100%', 'height': '5vh', 'align-items': 'center', 'margin-top': '1vh'})
    else:
        svg1 = None
        stats =  html.Div([
            html.P("#matches: " + str(len(siteLocator.matchedPeaks)), style = {'margin-left': '2vw'}),
            html.P("#shifted: " + str(len(siteLocator.shifted)), style = {'margin-left': '2vw'}),
        ], style = {'display': 'flex', 'flex-direction': 'row', 'flex-wrap': 'wrap', 'justify-content': 'left', 'width': '100%', 'height': '5vh', 'align-items': 'center', 'margin-top': '1vh'})
    
    svg2 = visualizer.highlightScores(mol1, scores)




    if (smiles2 is not None and len(smiles2) > 0):
        return dash_svg(svg1), dash_svg(svg2), "Success", base64.b64encode(pickle.dumps(siteLocator)).decode(), stats 
    return dash_svg(svg1), dash_svg(svg2), "Success", base64.b64encode(pickle.dumps(siteLocator)).decode(), stats

## update debugtext if click on bar chart
@app.callback(
    Output('peak_info', 'children'),
    Input('peaks', 'clickData'),
    State('siteLocatorObj', 'data'), prevent_initial_call=True)
def display_click_data(clickData, siteLocatorObj):
    if clickData:
        try:
            siteLocator = pickle.loads(base64.b64decode(siteLocatorObj))
            structures = siteLocator.get_structures_per_peak(float(clickData['points'][0]['x']))
            res = []
            for structure in structures:
                svg = dash_svg(visualizer.molToSVG(siteLocator.molMol, Chem.MolFromSmiles(structure, sanitize=False)))
                res.append(html.Img(src=svg, style={'width': '300px'}))
            return res
                
        except:
            return "siteLocator object not found"
    return None
   

@app.callback([ Output('USI1', 'value'),
                Output('USI2', 'value'),
                Output('SMILES1', 'value'),
                Output('SMILES2', 'value'),
                Output('options', 'value'),
                Output('filter_peaks_method', 'value'),
                Output('filter_peaks_variable', 'value'),
                Output('mz_tolerance', 'value'),
                Output('ppm', 'value')],
                Input('url', 'href'), prevent_initial_call=True)
def _content(href: str):
    print("this got called", href)
    try:
        f = furl(href)
        USI1 = f.args.get('USI1', 'mzspec:GNPS:TASK-5700dee92610412ea452a4262add2b93-f.MSV000086107/ccms_peak/VVP3-2_EtAc_MeOh.mzML:scan:2303')
        USI2 = f.args.get('USI2', 'mzspec:GNPS:TASK-5700dee92610412ea452a4262add2b93-f.MSV000086107/ccms_peak/VVP77_MeOh.mzML:scan:3024')
        SMILES1 = f.args.get('SMILES1', 'OC1=CC=CC=C1C2=NC(C3N(C)C(C(O)=O)CS3)CS2')
        SMILES2 = f.args.get('SMILES2', 'OC1=CC=CC=C1C2=NC(C3N(C)C(C(OC)=O)CS3)CS2')
        ## read arguments from url
        ppm = f.args.get('ppm', 1.01)
        mz_tolerance = f.args.get('mz_tolerance', 0.1)
        filter_peaks_method = f.args.get('filter_peaks_method', 'top_k')
        filter_peaks_variable = f.args.get('filter_peaks_variable', 50)

        options = f.args.get('options', [])
        if options is not None and type(options) is str:
            options = json.loads(options.replace("'", '"'))

        return USI1, USI2, SMILES1, SMILES2, options, filter_peaks_method, filter_peaks_variable, mz_tolerance, ppm
    except:
        print("error in url")
        return None, None, None, None, None, None, None, None, None

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

    app.run_server(debug = True)#, port=80, host="0.0.0.0")