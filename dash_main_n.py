from dash import Dash, html, dcc, Input, Output, State
from Dash_interface import chart_section_n, input_section_n, url_manager_n, computation_n
from arguments import args
import ModificationSiteLocator as modSite
import Compound_n as compound
import visualizer as vis
import utils_n as utils
import handle_network as hn
import dash_bootstrap_components as dbc
from flask import Flask


server = Flask(__name__)
app = Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])
app.title = 'Modification Site Locator'
server = app.server

if __name__ == '__main__':
    layout = [dcc.Store(id='siteLocatorObj'), dcc.Store(id='InputData'), dcc.Store(id='urlData'), html.P(id='ErrorMessage')]
    # layout.append(url_manager_n.get_layout())
    layout.append(input_section_n.get_layout(args, hn, utils))
    layout.append(chart_section_n.get_layout())
    app.layout = html.Div(layout)

    # url_manager_n.get_callbacks(app, args)
    input_section_n.get_callbacks(app)
    computation_n.get_callbacks(app, modSite, compound, hn)
    chart_section_n.get_callbacks(app, vis, utils)

    app.run_server(debug=False, port=5000, host="0.0.0.0") # app.run_server(debug = True, port=8081)#port=80, host="0.0.0.0")