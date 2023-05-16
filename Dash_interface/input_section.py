from dash import Dash, html, dcc, Input, Output, State, dash_table
import dash_daq as daq
import dash_bootstrap_components as dbc

def get_layout(passedArgs):
    print ("passedArgs", passedArgs)
    args = {}
    for key in passedArgs:
        try:
            args[key] = passedArgs[key].get('defult')
        except:
            args[key] = passedArgs[key]

    layer1 = []
    for item in ['USI1', 'USI2', 'SMILES1', 'SMILES2']:
        layer1.append(dbc.InputGroup(
            [dbc.InputGroupText(item), dbc.Input(placeholder=item,id=item, value = args.get(item, ""))],
            style = {'width': '45vw', 'margin': '1vh'}
        ))
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