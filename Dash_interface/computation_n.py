from dash import Dash, html, dcc, Input, Output, State, dash_table
import base64
import pickle

def get_callbacks(app, modSite, Compound, hn):
    
    @app.callback(
        [Output('siteLocatorObj', 'data'), Output('ErrorMessage', 'children')],
        Input('InputData', 'data'),
        )
    def calculate_module(data):
        if data is None:
            return None, None
        main_info = hn.getDataFromUsi(data['USI1'])
        main_info['Adduct'] = "M+H"
        mod_info = hn.getDataFromUsi(data['USI2'])
        mod_info['Adduct'] = "M+H"
        main_compound = Compound.Compound(main_info, data['SMILES1'])
        mod_compound = Compound.Compound(mod_info, data.get('SMILES2', None))
        siteLocator = modSite.ModificationSiteLocator(main_compound, mod_compound)
        if siteLocator.main_compound.Precursor_MZ > siteLocator.modified_compound.Precursor_MZ:
            return None, "Molecule precursor mass is higher than modified precursor mass"
        else:
            return base64.b64encode(pickle.dumps(siteLocator)).decode(), "Success"