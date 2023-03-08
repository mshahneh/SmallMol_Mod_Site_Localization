from dash import Dash, html, dcc, Input, Output, State, dash_table
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
print(columns)

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
        html.P(id='table_out', children=len(df)),
        # refresh button
        html.Button('Refresh', id='refresh', n_clicks=0),
        dash_table.DataTable(
        id='my_table',
        columns=columns,
        data=df.to_dict('records'),
        style_cell={
        'overflow': 'hidden',
        'textOverflow': 'ellipsis',
        'textAlign': 'left',
        },
        style_data = {        'whiteSpace': 'normal',
        'height': 'auto',},
        markdown_options={"html": True, "link_target": "_self"},
        tooltip_delay=0,
        tooltip_duration=None,
        tooltip_header={i['name']: i['name'] for i in columns},
        page_current=0,
        page_size=10,
        page_action='custom',
        row_selectable="single",
        sort_action='custom',
        sort_mode='single',
        sort_by=[]),
        html.Img(id = 'prediction'),
        html.Img(id = 'substr'),
        #bar chart
        dcc.Graph(id='peaks'),
        # text
        html.Div(id='debug_text'),
    ])

@app.callback(
    [
        Output('substr', 'src'),
        Output('prediction', 'src'),
        Output('peaks', 'figure')
    ],
    Input('my_table', 'selected_rows'),
    State('my_table', 'data'))
def update_graphs(active_rows, local_df):
    if active_rows:
        print ("clicked on row:", active_rows[0])
        m0 = local_df[active_rows[0]]["mol1ID"]
        m1 = local_df[active_rows[0]]["mol2ID"]
        siteLocator = modSite.SiteLocator(utils.generate_usi(m0, library), utils.generate_usi(m1, library), cachedStructures[m0])
        scores_unshifted, scores_shifted = siteLocator.calculate_score()
        scores = siteLocator.distance_score(scores_unshifted, scores_shifted)
        if (max(scores.values()) > 0):
            scores = [scores[x] / max(scores.values()) for x in scores.keys()]
        else:
            scores = [0 for x in scores.keys()]
        svg1 = visualizer.molToSVG(cachedStructures[m0], cachedStructures[m1], True)
        svg2 = visualizer.highlightScores(cachedStructures[m0], scores)

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

        print(max(y1), len(y1), min(y2), len(y2))

        # get index of top 10 peaks in y2
        topPeaksInx = sorted(range(len(y2)), key=lambda i: y2[i])[-topPeakCount:]
        x2 = [round(x2[i], 2) for i in topPeaksInx]
        y2 = [y2[i] for i in topPeaksInx]
        y2 = [-x / max(y2)*100 for x in y2]
        print(x1, y1, x2, y2)
        fig.add_trace(go.Bar(x=x2, y=y2, name='modified molecule peaks', width=2))
        fig.update_layout(
            title="Peaks",
            bargap=0,
            xaxis_title="m/z",
            yaxis_title="intensity",
            legend_title="Peak Type",
        )

        # save siteLocator object
        with(open("tempSiteLocator.pkl", "wb")) as f:
            pickle.dump(siteLocator, f)

        return dash_svg(svg1), dash_svg(svg2), fig
    return None, None, go.Figure()

## update debugtext if click on bar chart
@app.callback(
    Output('debug_text', 'children'),
    Input('peaks', 'clickData'))
def display_click_data(clickData):
    if clickData:
        try:
            with(open("tempSiteLocator.pkl", "rb")) as f:
                siteLocator = pickle.load(f)
            structures = siteLocator.get_structures_per_peak(float(clickData['points'][0]['x']))
            res = []
            for structure in structures:
                svg = dash_svg(visualizer.molToSVG(siteLocator.molMol, Chem.MolFromSmiles(structure)))
                res.append(html.Img(src=svg, style={'width': '300px'}))
            return res
                
        except:
            return "siteLocator object not found"
    return None


@app.callback(
    Output('my_table', 'data'),
    Input('my_table', "page_current"),
    Input('my_table', "page_size"),
    Input('my_table', 'sort_by'))
def update_table(page_current, page_size, sort_by):
    if len(sort_by):
        dff = df.sort_values(
            sort_by[0]['column_id'],
            ascending=sort_by[0]['direction'] == 'asc',
            inplace=False
        )
    else:
        # No sort is applied
        dff = df

    temp = dff.iloc[page_current*page_size:(page_current+ 1)*page_size].to_dict('records')
    for i in range(len(temp)):
        svg = visualizer.molToSVG(cachedStructures[temp[i]['mol1ID']], cachedStructures[temp[i]['mol2ID']], True)
        temp[i]['structure'] = f"<img src=\"{dash_svg(svg)}\" style=\"height:100px;width:100px;\">"

    return temp
        

def dash_svg(text):
    """
    Generates a svg link for dash from a avg text.
    """
    svg = base64.b64encode(text.encode('utf-8'))
    return 'data:image/svg+xml;base64,{}'.format(svg.decode())

if __name__ == '__main__':
    #read data

    for match in tqdm(matches[1]):
        # break if we have enough matches
        if len(df) > 15:
            break
        try:
            m0, m1 = match
            molMol = cachedStructures[m1]
            modifMol = cachedStructures[m0]
            molUsi = utils.generate_usi(m1, data_dict_filtered[m1]['library_membership'])
            modifUsi = utils.generate_usi(m0, data_dict_filtered[m0]['library_membership'])
            molSmiles = data_dict_filtered[m1]['Smiles']
            site = modSite.SiteLocator(molUsi, modifUsi, molSmiles)
            modifLoc = list(utils.calculateModificationSites(modifMol, molMol, False))
            df = pd.concat([df, pd.DataFrame.from_records([{"mol1ID": m0, "mol2ID": m1, "score": site.accuracy_score(modifLoc[0]), 
            "# matched peaks": len(site.matchedPeaks), "# shifted peaks": len(site.shifted), "# unshifted peaks": len(site.unshifted)}])])
        except:
            traceback.print_exc()
    
    app.run_server()#port=80, host="0.0.0.0")