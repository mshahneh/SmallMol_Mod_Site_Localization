from dash import Dash, html, dcc, Input, Output, State, dash_table
import dash_daq as daq
import dash_bootstrap_components as dbc
import pandas as pd
from dash.dash_table.Format import Format, Scheme, Trim
import pickle
import traceback
from tqdm import tqdm

def get_layout(passedArgs, hn, utils):

    myData = [{'USI1': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010126753', 'USI2': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010118619', 'Smiles1': 'O=C(O)[C@@H]1CCC(O)=N1', 'Smiles2': 'CCOC(=O)C1CCC(O)=N1'}, {'USI1': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010125637', 'USI2': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010106370', 'Smiles1': 'COc1cc(C(=O)O)ccc1O', 'Smiles2': 'CCOC(=O)c1ccc(O)c(OC)c1'}, {'USI1': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010125095', 'USI2': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010120285', 'Smiles1': 'O=C(O)c1cc(O)c(O)c(O)c1', 'Smiles2': 'CCOC(=O)c1cc(O)c(O)c(O)c1'}, {'USI1': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010122676', 'USI2': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010108770', 'Smiles1': 'COc1cc(C=O)ccc1O', 'Smiles2': 'COc1cc(C=O)ccc1OC(C)C'}, {'USI1': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010122675', 'USI2': 'mzspec:GNPS:BERKELEY-LAB:accession:CCMSLIB00010104696', 'Smiles1': 'COc1cc(C=O)ccc1O', 'Smiles2': 'CCCOc1ccc(C=O)cc1OC'}]
    columns=["USI1", "Smiles1", "USI2", "Smiles2"]
    df = pd.DataFrame(myData, columns=["USI1", "Smiles1", "USI2", "Smiles2"])

    print ("passedArgs", passedArgs)
    args = {}
    for key in passedArgs:
        try:
            args[key] = passedArgs[key].get('defult')
        except:
            args[key] = passedArgs[key]
    
    # table layer
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
    options = [{'label': 'Presence Only', 'value': 'presence_only'}, {'label': 'Consider Intensity', 'value': 'consider_intensity'}]#, {'label': 'Shifted Only', 'value': 'shifted_only'}]
    for item in options:
        myOptions.append(daq.BooleanSwitch(
            id=item['value'],
            on=args.get(item['value'], False),
            label=item['label'],
            labelPosition="Left" , style={'display': 'inline-block', 'margin': '2vh 1vw'})
        )

    optionArguments = []
    # optionArguments = [dcc.Dropdown(['intensity', 'top_k', 'none'], 'top_k', id='filter_peaks_method')]
    for item in [ 'mz_tolerance']:#, 'filter_peaks_variable', 'ppm']:
        optionArguments.append(dbc.InputGroup(
            [dbc.InputGroupText(item), dbc.Input(placeholder=item,id=item, value = args.get(item, ""))],
            style={"height": "2vh"}
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
    with open('data/libraries/BERKELEY-LAB/cachedStructures.pkl', 'rb') as f:
        cachedStructures = pickle.load(f)
    @app.callback(
        Output('inputs', 'children'),
        Input('urlData', 'data'),
        prevent_initial_call=True
    )
    def _content(data):
        print(data)
        return get_layout(data)

    @app.callback(
        Output('InputData', 'data'),
        Input('update', 'n_clicks'),
        State('InputData', 'data'),
        State('USI1', 'value'),
        State('USI2', 'value'),
        State('SMILES1', 'value'),
        State('SMILES2', 'value'),
        # State('filter_peaks_method', 'value'),
        # State('filter_peaks_variable', 'value'),
        State('mz_tolerance', 'value'),
        # State('ppm', 'value'),
        State('presence_only', 'on'),
        State('consider_intensity', 'on'),
        # State('shifted_only', 'on'),
        prevent_initial_call=True
    )
    def _content(n_clicks, data, USI1, USI2, SMILES1, SMILES2, mz_tolerance, presence_only, consider_intensity):
        temp =  {
            'USI1': USI1,
            'USI2': USI2,
            'SMILES1': SMILES1,
            'SMILES2': SMILES2,
            # 'filter_peaks_method': filter_peaks_method,
            # 'filter_peaks_variable': filter_peaks_variable,
            'mz_tolerance': mz_tolerance,
            # 'ppm': ppm,
            'presence_only': presence_only,
            'consider_intensity': consider_intensity,
            # 'shifted_only': shifted_only,
        }
        if data is None:
            data = {}
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