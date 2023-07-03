from dash import Dash, html, dcc, Input, Output, State
from Dash_interface import chart_section, input_section, url_manager, computation
from arguments import args
import SiteLocator as modSite
import visualizer as vis
import utils as utils
import handle_network as hn
import dash_bootstrap_components as dbc

app = Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

if __name__ == '__main__':
    layout = [dcc.Store(id='siteLocatorObj'), dcc.Store(id='InputData'), dcc.Store(id='urlData'), html.P(id='ErrorMessage')]
    layout.append(url_manager.get_layout())
    layout.append(input_section.get_layout(args, hn, utils))
    layout.append(chart_section.get_layout())
    app.layout = html.Div(layout)

    url_manager.get_callbacks(app, args)
    input_section.get_callbacks(app)
    computation.get_callbacks(app, modSite)
    chart_section.get_callbacks(app, vis, utils)

    app.run_server(debug = True, port=80, host="0.0.0.0")