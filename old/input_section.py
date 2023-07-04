from dash import Dash, html, dcc, Input, Output, State, dash_table
import dash_daq as daq
import dash_bootstrap_components as dbc
import pandas as pd
from dash.dash_table.Format import Format, Scheme, Trim
import pickle
import traceback
from tqdm import tqdm

def get_layout(passedArgs, hn, utils):
    # df = pd.DataFrame(columns=["USI1", "Smiles1", "USI2", "Smiles2"])
    # columns=[{"name": i, "id": i} for i in df.columns]
    # # columns.append({"name": "structure", "id": "structure", "presentation": "markdown"})
    # # columns[2]["format"] = Format(precision=2, scheme=Scheme.fixed)

    # library = "BERKELEY-LAB"
    # with open('data/libraries/{}/data_dict_filtered.pkl'.format(library), 'rb') as f:
    #     data_dict_filtered = pickle.load(f)

    # with open('data/libraries/{}/matches.pkl'.format(library), 'rb') as f:
    #     matches = pickle.load(f)

    # with open('data/libraries/{}/cachedStructures.pkl'.format(library), 'rb') as f:
    #     cachedStructures = pickle.load(f)
    # for match in tqdm(matches[1]):
    #     # break if we have enough matches
    #     if len(df) > 15:
    #         break
    #     try:
    #         m0, m1 = match
    #         molMol = cachedStructures[m1]
    #         modifMol = cachedStructures[m0]
    #         molUsi = hn.generate_usi(m1, data_dict_filtered[m1]['library_membership'])
    #         modifUsi = hn.generate_usi(m0, data_dict_filtered[m0]['library_membership'])
    #         molSmiles = data_dict_filtered[m1]['Smiles']
    #         modifSmiles = data_dict_filtered[m0]['Smiles']
    #         modifLoc = utils.calculateModificationSites(modifMol, molMol, False)
    #         df = pd.concat([df, pd.DataFrame.from_records([{"USI1": molUsi, "Smiles1": molSmiles, "USI2": modifUsi, "Smiles2": modifSmiles}])])
    #     except:
    #         traceback.print_exc()
    
    # # select 5 rows from the dataframe
    # df = df.iloc[:5]

    myData = [{'USI1': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010126753', 'USI2': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010118619', 'Smiles1': 'O=C(O)[C@@H]1CCC(O)=N1', 'Smiles2': 'CCOC(=O)C1CCC(O)=N1'}, {'USI1': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010125637', 'USI2': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010106370', 'Smiles1': 'COc1cc(C(=O)O)ccc1O', 'Smiles2': 'CCOC(=O)c1ccc(O)c(OC)c1'}, {'USI1': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010125095', 'USI2': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010120285', 'Smiles1': 'O=C(O)c1cc(O)c(O)c(O)c1', 'Smiles2': 'CCOC(=O)c1cc(O)c(O)c(O)c1'}, {'USI1': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010122676', 'USI2': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010108770', 'Smiles1': 'COc1cc(C=O)ccc1O', 'Smiles2': 'COc1cc(C=O)ccc1OC(C)C'}, {'USI1': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010122675', 'USI2': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010104696', 'Smiles1': 'COc1cc(C=O)ccc1O', 'Smiles2': 'CCCOc1ccc(C=O)cc1OC'}]
    df = pd.DataFrame(myData, columns=["USI1", "Smiles1", "USI2", "Smiles2"])
    columns=[{"name": i, "id": i} for i in df.columns]
    print ("passedArgs", passedArgs)
    args = {}
    for key in passedArgs:
        try:
            args[key] = passedArgs[key].get('defult')
        except:
            args[key] = passedArgs[key]
    
    table =         dash_table.DataTable(
        id='my_table',
        columns=columns,
        data=df.to_dict('records'),
        style_cell={
        'overflow': 'hidden',
        'textOverflow': 'ellipsis',
        'textAlign': 'left',
        'maxWidth': '25%',
        },
        style_data = {        'whiteSpace': 'normal',
        'height': 'auto',},
        markdown_options={"html": True, "link_target": "_self"},
        tooltip_delay=0,
        tooltip_duration=None,
        tooltip_header={i['name']: i['name'] for i in columns},
        page_current=0,
        row_selectable="single",
        style_table={ 'width': '95%', 'margin': 'auto'},
        )
    inp1 = [html.H5('Base Compound',style = {'width': '100%', 'margin': '1vh'})]
    for item in ['USI1', 'SMILES1']:
        inp1.append(dbc.InputGroup(
            [dbc.InputGroupText(item), dbc.Input(placeholder=item,id=item, value = args.get(item, ""))],
            style = {'width': '45vw', 'margin': '1vh'}
        ))
    inp2 = [html.H5('Modified Compound',style = {'width': '100%', 'margin': '1vh'})]
    for item in ['USI2', 'SMILES2']:
        inp2.append(dbc.InputGroup(
            [dbc.InputGroupText(item), dbc.Input(placeholder=item,id=item, value = args.get(item, ""))],
            style = {'width': '45vw', 'margin': '1vh'}
        ))
    layer1 = []
    layer1.append(html.Div(children = inp1, style = {'border': '1px solid #aabbdd', 'display': 'flex', 'flex-direction': 'column', 'justify-content': 'center', 'align-items': 'center', 'padding':'10px'}))
    layer1.append(html.Div(children = inp2, style = {'border': '1px solid #aabbdd', 'display': 'flex', 'flex-direction': 'column', 'justify-content': 'center', 'align-items': 'center', 'padding':'10px'}))
    layer1 = html.Div(children = layer1, style = {'display': 'flex', 'flex-direction': 'row', 'flex-wrap': 'wrap', 'justify-content': 'space-around', 'width': '100%', 'min-height': '5vh'})

    myOptions = []
    for item in ['presence_only', 'consider_intensity', 'shifted_only']:
        myOptions.append(daq.BooleanSwitch(
            id=item,
            on=args.get(item, False),
            label=item,
            labelPosition="top", style={'display': 'inline-block', 'margin-right': '1vw'})
        )

    optionArguments = [dcc.Dropdown(['intensity', 'top_k', 'none'], 'top_k', id='filter_peaks_method')]
    for item in ['filter_peaks_variable', 'mz_tolerance', 'ppm']:
        optionArguments.append(dbc.InputGroup(
            [dbc.InputGroupText(item), dbc.Input(placeholder=item,id=item, value = args.get(item, ""))],
            style={'width':'30%',"height": "2vh"}
        ))

    return html.Div(id = 'inputs', children = 
        [
            html.Div(children=[table], style = {'width': '100%', 'margin-bottom': '1vh'}),
            layer1,
            
            html.Div(children = [
                html.Div(
                        id='options',
                        children=myOptions,
                    style={'display': 'flex', 'flex-direction': 'row','margin-right': '1vw'}),
                html.Div(id = 'arguments', children = optionArguments, style = {'display': 'flex', 'flex-direction': 'row', 'justify-content': 'center', 'min-height': '5vh', 'align-items': 'center', 'margin-top': '1vh'}),
            ], style = {'display': 'flex', 'flex-direction': 'row', 'justify-content': 'center', 'align-items': 'center'}),
            
            # html.Div(children = [
            #     dcc.Input(
            #         id="helperUSI",
            #         type="text",
            #         placeholder="helperUSI",
            #         style={'width':'23%'}
            #     ),
            #     dcc.Input(
            #         id="helperSMILES",
            #         type="text",
            #         placeholder="helperSMILES",
            #         style={'width':'23%'}
            #     ),
            #         # html.Button('update', id='updateHelper', n_clicks=0)
            # ], style={'display': 'hidden'}
            # ),
            dbc.Button('update', id='update', n_clicks=0)
        ], style = {'display': 'flex', 'flex-direction': 'column', 'justify-content': 'center', 'align-items': 'center', 'width': '100%','margin-bottom': '1vh'})

def get_callbacks(app):
    # @app.callback(
    #     Output('inputs', 'children'),
    #     Input('urlData', 'data'),
    #     prevent_initial_call=True
    # )
    # def _content(data):
    #     print(data)
    #     return get_layout(data)

    @app.callback(
        Output('InputData', 'data'),
        Input('update', 'n_clicks'),
        State('InputData', 'data'),
        State('USI1', 'value'),
        State('USI2', 'value'),
        State('SMILES1', 'value'),
        State('SMILES2', 'value'),
        State('filter_peaks_method', 'value'),
        State('filter_peaks_variable', 'value'),
        State('mz_tolerance', 'value'),
        State('ppm', 'value'),
        State('presence_only', 'on'),
        State('consider_intensity', 'on'),
        State('shifted_only', 'on'),
        prevent_initial_call=True
    )
    def _content(n_clicks, data, USI1, USI2, SMILES1, SMILES2, filter_peaks_method, filter_peaks_variable, mz_tolerance, ppm, presence_only, consider_intensity, shifted_only):
        temp =  {
            'USI1': USI1,
            'USI2': USI2,
            'SMILES1': SMILES1,
            'SMILES2': SMILES2,
            'filter_peaks_method': filter_peaks_method,
            'filter_peaks_variable': filter_peaks_variable,
            'mz_tolerance': mz_tolerance,
            'ppm': ppm,
            'presence_only': presence_only,
            'consider_intensity': consider_intensity,
            'shifted_only': shifted_only,
        }
        data.update(temp)
        return data

    @app.callback(
    [
        Output('USI1', 'value'),
        Output('USI2', 'value'),
        Output('SMILES1', 'value'),
        Output('SMILES2', 'value'),
    ],
    Input('my_table', 'selected_rows'),
    State('my_table', 'data'))
    def update_inputs_from_table(active_rows, local_df):
        if active_rows is None or len(active_rows) == 0:
            return "", "", "", ""
        row = local_df[active_rows[0]]
        return row['USI1'], row['USI2'], row['Smiles1'], row['Smiles2']