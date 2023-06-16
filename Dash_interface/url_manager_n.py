from dash import Dash, html, dcc, Input, Output, State, dash_table
from furl import furl
from copy import deepcopy
import dash_bootstrap_components as dbc

def get_layout():
    return html.Div([dcc.Location(id='url', refresh=False),
                      dbc.Button('share', id='share',outline=True, color="secondary", n_clicks=0, style={"position":"absolute", "right":"10px", "top":"10px"}),
                      html.Div(id="share_dialog", children=[html.P(id='share_url'),
                                                             html.Button('x', id='close_share_dialog', style={"position":"absolute", "right":"10px", "top":"10px"})
                                                            ], style={"display":"none"})])

def get_callbacks(app, set_up_arguments):
    @app.callback([Output('urlData', 'data'),
                   Output('InputData', 'data', allow_duplicate=True),
                    Output('ErrorMessage', 'children', allow_duplicate=True)],
                    Input('url', 'href'), prevent_initial_call=True)
    def _content(href: str):
        try:
            f = furl(href)
            data = {}
            for argument in set_up_arguments:
                data[set_up_arguments[argument]['name']] = f.args.get(set_up_arguments[argument]['short'], set_up_arguments[argument]['defult'])
            print (data)
            # copy data to InputData
            return data, deepcopy(data), "Success"
        except:
            return {}, "Error in URL"

    @app.callback([Output('share_url', 'children'),
                    Output('share_dialog', 'style',  allow_duplicate=True)],
                    Input('share', 'n_clicks'),
                    [State('InputData', 'data'),
                     State('url', 'href')], prevent_initial_call=True)
    def _share(n_clicks, data, href):
        if n_clicks == 0:
            return "Share URL", {'display': 'none'}
        else:
            f = furl(href)
            # clear all arguments
            f.args = {}
            for argument in set_up_arguments:
                if data[set_up_arguments[argument]['name']] != set_up_arguments[argument]['defult']:
                    f.args[set_up_arguments[argument]['short']] = data[set_up_arguments[argument]['name']]
            return f.url, {"position":"absolute", 'display': 'block', "width":"80vw", "height":"20vh", "top":"10vh", "left":"10vw", 
                           "background-color":"rgba(0,0,0,0.8)", "z-index":"3", "padding":"15px", "color":"white"}
    
    @app.callback(Output('share_dialog', 'style'),
                    Input('close_share_dialog', 'n_clicks'), prevent_initial_call=True)
    def _close_share_dialog(n_clicks):
        return {'display': 'none'}
    
    