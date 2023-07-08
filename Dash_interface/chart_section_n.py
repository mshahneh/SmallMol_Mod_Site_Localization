from dash import Dash, html, dcc, Input, Output, State, dash_table
import pickle
import base64
import plotly.graph_objects as go
import rdkit.Chem as Chem
import dash_bootstrap_components as dbc
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
            html.Div(id = 'stats'),
            html.Div(id = "draw_molecule_result", style = {'display': 'flex', 'flex-direction': 'row', 'justify-content': 'center', 'align-items': 'center'})
        ],  style = {'display': 'flex', 'flex-direction': 'row', 'justify-content': 'space-around', 'align-items': 'center', 'position': 'relative', 'z-index': '2'}),

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
        [Output ('draw_molecule_result', 'children'),
        Output('stats', 'children')],
        Input('siteLocatorObj', 'data'),
        State('InputData', 'data')
    )
    def update_stats(siteLocatorObj, args):
        if (siteLocatorObj == None):
            return None, None
        siteLocator = pickle.loads(base64.b64decode(siteLocatorObj))
        scores = siteLocator.generate_probabilities()
        stats = []
        stats.append(("number of matched peaks" , str(len(siteLocator.matched_peaks))))
        stats.append(("number of shifted" , str(len(siteLocator.shifted))))
        stats.append(("delta weight", str(round(abs(siteLocator.main_compound.Precursor_MZ - siteLocator.modified_compound.Precursor_MZ), 4))))
        draw_output = []
        if (siteLocator.modified_compound.structure is not None):
            svg1 = visualizer.molToSVG(siteLocator.modified_compound.structure, siteLocator.main_compound.structure, True)
            trueSite = utils.calculateModificationSites(siteLocator.modified_compound.structure, siteLocator.main_compound.structure, False)[0]
            stats.append(("Score:", str(round(siteLocator.calculate_score(trueSite, "temp" , scores), 4))))
            draw_output.append(html.Div(id = 'prediction_section', children = [
                    html.H3("Modification Highlight "),
                     html.Img(src = dash_svg(svg1),
                              style={'max-width': '30vw'})
            ], style={'margin-right': '1vw', 'padding-right': '1vw', 'border-right': '1px dashed #aaa'}))

        svg2 = visualizer.highlightScores(siteLocator.main_compound.structure, scores)
        draw_output.append(html.Div(id = 'substr_section', children = [
                    html.H3("Prediction"),
                    html.Img(src = dash_svg(svg2),
                             style={'max-width': '30vw'})
            ]))


        stats_section = []
        for stat in stats:
            print("stat", stat)
            stats_section.append(html.Div(children = [
                dbc.Badge(stat[0], color="info", className="me-1"),
                html.P(children = stat[1], style={'margin': '0px', 'margin-left': '5px'}),
            ], style = {'display': 'flex', 'flex-direction': 'row','align-items': 'center', 'justify-content': 'space-between', 'background-color': '#f8f9fa', 'margin': '3px', 'padding': '0px', 'border-radius': '5px'}))
        
        stats = dbc.Card(children = stats_section, style = {'flex':'1', 'width': '100%', 'padding': '10px'})

        
        return draw_output, stats

    
    @app.callback(
        [Output('peaks', 'figure'),
        Output('peaks-component', 'style')],
        Input('siteLocatorObj', 'data'),
        Input('peak_filter_slider', 'value'), prevent_initial_call=True)
    def update_peaks(data, slider_value):
        print("update_peaks")
        if (data == None):
            return {}, {'display': 'none'}
        siteLocator = pickle.loads(base64.b64decode(data))
        fig = go.Figure()
        typesInxMain = {'matched_shifted': [], 'matched_unshifted': [], 'unmatched': []}
        x1 = []
        y1 = []
        for peak in siteLocator.main_compound.peaks:
            x1.append(peak[0])
            y1.append(peak[1])
        

        topPeakCount = slider_value
        topPeaksInxModif = sorted(range(len(y1)), key=lambda i: y1[i])[-topPeakCount:]
        for i in topPeaksInxModif:
            flag = False
            for j in siteLocator.matched_peaks:
                if (j[0] == i):
                    if (abs(siteLocator.main_compound.peaks[i][0] - siteLocator.modified_compound.peaks[j[1]][0]) > siteLocator.args['mz_tolerance']):
                        typesInxMain['matched_shifted'].append(i)
                    else:
                        typesInxMain['matched_unshifted'].append(i)
                    flag = True
                    break
            if (not flag):
                typesInxMain['unmatched'].append(i)

        typesInxModified = {'matched_shifted': [], 'matched_unshifted': [], 'unmatched': []}
        x2 = []
        y2 = []
        for peak in siteLocator.modified_compound.peaks:
            x2.append(peak[0])
            y2.append(peak[1])
        
        topPeaksInxModif = sorted(range(len(y2)), key=lambda i: y2[i])[-topPeakCount:]
        for i in topPeaksInxModif:
            flag = False
            for j in siteLocator.matched_peaks:
                if (j[1] == i):
                    if (abs(siteLocator.main_compound.peaks[j[0]][0] - siteLocator.modified_compound.peaks[j[1]][0]) > 0.1):
                        typesInxModified['matched_shifted'].append([i, j[0]])
                    else:
                        typesInxModified['matched_unshifted'].append([i, j[0]])
                    flag = True
                    break
            if (not flag):
                typesInxModified['unmatched'].append([i, -1])
        
        for inx_type in typesInxMain:
            x_main = [round(x1[j], 4) for j in typesInxMain[inx_type]]
            y1_ = [y1[j] for j in typesInxMain[inx_type]]
            y_main = [x / max(y1_) * 100 for x in y1_]
            x_modified = [round(x2[j[0]], 4) for j in typesInxModified[inx_type]]
            y2_ = [y2[j[0]] for j in typesInxModified[inx_type]]
            y_modified = [-j / max(y2_) * 100 for j in y2_]
            indicis = typesInxMain[inx_type]
            x_ = x_main + x_modified
            y_ = y_main + y_modified
            if (inx_type == 'unmatched'):
                fig.add_trace(go.Bar(x=x_, y=y_, hovertext=indicis, name=inx_type, width=0.2, visible='legendonly'))
            else:
                fig.add_trace(go.Bar(x=x_, y=y_, hovertext=indicis, name=inx_type, width=0.2))
        
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
                clicked_peak_x = float(clickData['points'][0]['x'])
                clicked_peak_y = float(clickData['points'][0]['y'])
                if (clicked_peak_y < 0):
                    return "only base compound peaks are clickable"
                structures, result_posibility_indicies = siteLocator.get_structures_by_peak_weight(float(clickData['points'][0]['x']))
                res = []
                for posibility_index in result_posibility_indicies:
                    svg = dash_svg(visualizer.highlightMolIndices(siteLocator.main_compound.structure, posibility_index))
                    res.append(html.Img(src=svg, style={'width': '300px'}))
                
                return res
            except:
                import traceback
                traceback.print_exc()
                return "siteLocator object not found"
        return None
