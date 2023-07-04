from dash import Dash, html, dcc, Input, Output, State, dash_table
import pickle
import base64
import plotly.graph_objects as go
import rdkit.Chem as Chem
import numpy as np

def dash_svg(text):
    """
    Generates a svg link for dash from a avg text.
    """
    try:
        svg = base64.b64encode(text.encode('utf-8'))
        return 'data:image/svg+xml;base64,{}'.format(svg.decode())
    except:
        return None

def get_layout():
     return html.Div(id = "results", children = [
        html.Div(children = [
            html.Div(id = 'stats', style={"margin-right": "1vw"}),
            html.Div(children = [
                html.Img(id = 'prediction', style = {'width': '50%'}),
                html.Img(id = 'substr', style = {'width': '50%'}),
            ], style = {'display': 'flex', 'flex-direction': 'row', 'justify-content': 'center', 'align-items': 'center'})
        ],  style = {'display': 'flex', 'flex-direction': 'row', 'justify-content': 'center', 'align-items': 'center', 'position': 'relative', 'z-index': '2'}),

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
        html.Div(id='peak_info')
    ])

def get_callbacks(app, visualizer, utils):
    @app.callback(
        [Output('prediction', 'src'),
        Output('substr', 'src'),
        Output('stats', 'children')],
        Input('siteLocatorObj', 'data'),
        State('InputData', 'data')
    )
    def update_stats(siteLocatorObj, args):
        if (siteLocatorObj == None):
            return None, None, None
        siteLocator = pickle.loads(base64.b64decode(siteLocatorObj))
        scores_unshifted, scores_shifted = siteLocator.calculate_score(peak_presence_only = args['presence_only'], consider_intensity = args['presence_only'])
        print ("debugging: scores shifted:", scores_unshifted, scores_shifted, (args['shifted_only']==False), args['presence_only'], args['presence_only'])
        scores = siteLocator.distance_score(scores_unshifted, scores_shifted, combine = (args['shifted_only']==False))
    
        isMax = None
        score = None

        smiles1 = args['SMILES1']
        smiles2 = args['SMILES2']
        mol1 = Chem.MolFromSmiles(smiles1)

        if (smiles2 is not None and len(smiles2) > 0):
            svg1 = visualizer.molToSVG(Chem.MolFromSmiles(smiles2), mol1, True)
            modifLoc = list(utils.calculateModificationSites(Chem.MolFromSmiles(smiles2), mol1, False))
            accuracy_score = siteLocator.accuracy_score(modifLoc[0], peak_presence_only=args['presence_only'], combine=(args['shifted_only']==False), return_all=True, consider_intensity=args['consider_intensity'])
            isMax = accuracy_score['isMax']
            score = accuracy_score['score']
            stats =  html.Div([
                html.P("cosine: " + str(round(siteLocator.cosine, 4)), style = {'margin-left': '1.5vw'}),
                html.P("Score: " + str(round(score, 4)), style = {'margin-left': '1.5vw'}),
                html.P("is Max: " + str(isMax), style = {'margin-left': '1.5vw'}),
                html.P("#matches: " + str(len(siteLocator.matchedPeaks)), style = {'margin-left': '1.5vw'}),
                html.P("#shifted: " + str(len(siteLocator.shifted)), style = {'margin-left': '1.5vw'}),
                html.P("delta w:" + str(round(abs(siteLocator.molPrecursorMz - siteLocator.modifPrecursorMz), 4)), style = {'margin-left': '1.5vw'}),
            ], style = {'display': 'flex', 'flex-direction': 'row', 'flex-wrap': 'wrap', 'justify-content': 'left', 'width': '100%', 'height': '5vh', 'align-items': 'center', 'margin-top': '1vh'})
        else:
            svg1 = None
            stats =  html.Div([
                html.P("#matches: " + str(len(siteLocator.matchedPeaks)), style = {'margin-left': '2vw'}),
                html.P("#shifted: " + str(len(siteLocator.shifted)), style = {'margin-left': '2vw'}),
            ], style = {'display': 'flex', 'flex-direction': 'row', 'flex-wrap': 'wrap', 'justify-content': 'left', 'width': '100%', 'height': '5vh', 'align-items': 'center', 'margin-top': '1vh'})
        
        print("scores are:", scores)
        svg2 = visualizer.highlightScores(mol1, scores)
        return dash_svg(svg1), dash_svg(svg2), stats

    
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
                    if (abs(siteLocator.molPeaks[i][0] - siteLocator.modifPeaks[j[1]][0]) > siteLocator.args['mz_tolerance']):
                        typesInx['matched_shifted'].append(i)
                    else:
                        typesInx['matched_unshifted'].append(i)
                    flag = True
                    break
            if (not flag):
                typesInx['unmatched'].append(i)
        
        for i in typesInx:
            x1_ = [round(x1[j], 4) for j in typesInx[i]]
            y1_ = [y1[j] for j in typesInx[i]]
            y1_ = [x / max(y1_) * 100 for x in y1_]
            indicis = typesInx[i]
            if (i == 'unmatched'):
                fig.add_trace(go.Bar(x=x1_, y=y1_, hovertext=indicis, name=i, width=0.4, visible='legendonly'))
            else:
                fig.add_trace(go.Bar(x=x1_, y=y1_, hovertext=indicis, name=i, width=0.4))

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
                        typesInx['matched_shifted'].append([i, j[0]])
                    else:
                        typesInx['matched_unshifted'].append([i, j[0]])
                    flag = True
                    break
            if (not flag):
                typesInx['unmatched'].append([i, -1])
        

        for i in typesInx:
            x1_ = [round(x2[j[0]], 4) for j in typesInx[i]]
            y2_ = [y2[j[0]] for j in typesInx[i]]
            y2_ = [-j / max(y2_) * 100 for j in y2_]
            indicis = typesInx[i]
            if (i == 'unmatched'):
                fig.add_trace(go.Bar(x=x1_, y=y2_, hovertext=indicis, name=i, width=0.4, visible='legendonly'))
            else:
                fig.add_trace(go.Bar(x=x1_, y=y2_, hovertext=indicis, name=i, width=0.4))
        
        fig.update_layout(
            title="Peaks",
            bargap=0,
            xaxis_title="m/z",
            yaxis_title="intensity",
            legend_title="Peak Type",
        )
        return fig, {'display': 'flex', 'position': 'relative', 'flex-direction': 'column', 'justify-content': 'center', 'width': '100%', 'align-items': 'center', 'margin-top': '-4vh', "z-index": "1"}

    @app.callback(
        Output('peak_info', 'children'),
        Input('peaks', 'clickData'),
        State('siteLocatorObj', 'data'), prevent_initial_call=True)
    def display_click_data(clickData, siteLocatorObj):
        if clickData:
            try:
                siteLocator = pickle.loads(base64.b64decode(siteLocatorObj))
                structures, result_posibility_indicies = siteLocator.get_structures_per_peak(float(clickData['points'][0]['x']))
                res = []
                # for structure in structures:
                #     svg = dash_svg(visualizer.molToSVG(siteLocator.molMol, Chem.MolFromSmiles(structure, sanitize=False)))
                #     res.append(html.Img(src=svg, style={'width': '300px'}))
                
                print ("debugging in display click data", result_posibility_indicies)
                for posibility_index in result_posibility_indicies:
                    svg = dash_svg(visualizer.highlightMolIndices(siteLocator.molMol, posibility_index))
                    res.append(html.Img(src=svg, style={'width': '300px'}))
                
                return res
            except:
                return "siteLocator object not found"
        return None
