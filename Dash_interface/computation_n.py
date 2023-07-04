from dash import Dash, html, dcc, Input, Output, State, dash_table
import base64
import pickle

def get_callbacks(app, modSite, Compound, hn):
    
    @app.callback(
        [Output('siteLocatorObj', 'data'), Output('ErrorMessage', 'children')],
        Input('InputData', 'data'),
        )
    def calculate_module(data):
        print ('\n\n\n', data, '\n\n\n')
        if data is None:
            return None, None
        print("we actually get to here")
        main_info = hn.getDataFromUsi(data['USI1'])
        main_info['Adduct'] = "M+H"
        mod_info = hn.getDataFromUsi(data['USI2'])
        mod_info['Adduct'] = "M+H"
        main_compound = Compound.Compound(main_info, data['SMILES1'])
        mod_compound = Compound.Compound(mod_info, data['SMILES2'])
        siteLocator = modSite.ModificationSiteLocator(main_compound, mod_compound)
        if siteLocator.main_compound.Precursor_MZ > siteLocator.modified_compound.Precursor_MZ:
            return None, "Molecule precursor mass is higher than modified precursor mass"
        else:
            return base64.b64encode(pickle.dumps(siteLocator)).decode(), "Success"
    # scores_unshifted, scores_shifted = siteLocator.calculate_score(peak_presence_only = presense, consider_intensity = consider_intensity)
    # scores = siteLocator.distance_score(scores_unshifted, scores_shifted, combine = combine)
    
    # isMax = None
    # score = None
    # if (smiles2 is not None and len(smiles2) > 0):
    #     svg1 = visualizer.molToSVG(Chem.MolFromSmiles(smiles2), mol1, True)
    #     modifLoc = list(utils.calculateModificationSites(Chem.MolFromSmiles(smiles2), mol1, False))
    #     accuracy_score = siteLocator.accuracy_score(modifLoc[0], peak_presence_only=presense, combine=combine, return_all=True, consider_intensity=consider_intensity)
    #     isMax = accuracy_score['isMax']
    #     score = accuracy_score['score']
    #     stats =  html.Div([
    #         html.P("cosine: " + str(round(siteLocator.cosine, 4)), style = {'margin-left': '1.5vw'}),
    #         html.P("Score: " + str(round(score, 4)), style = {'margin-left': '1.5vw'}),
    #         html.P("is Max: " + str(isMax), style = {'margin-left': '1.5vw'}),
    #         html.P("#matches: " + str(len(siteLocator.matchedPeaks)), style = {'margin-left': '1.5vw'}),
    #         html.P("#shifted: " + str(len(siteLocator.shifted)), style = {'margin-left': '1.5vw'}),
    #         html.P("delta w:" + str(round(abs(siteLocator.molPrecursorMz - siteLocator.modifPrecursorMz), 4)), style = {'margin-left': '1.5vw'}),
    #     ], style = {'display': 'flex', 'flex-direction': 'row', 'flex-wrap': 'wrap', 'justify-content': 'left', 'width': '100%', 'height': '5vh', 'align-items': 'center', 'margin-top': '1vh'})
    # else:
    #     svg1 = None
    #     stats =  html.Div([
    #         html.P("#matches: " + str(len(siteLocator.matchedPeaks)), style = {'margin-left': '2vw'}),
    #         html.P("#shifted: " + str(len(siteLocator.shifted)), style = {'margin-left': '2vw'}),
    #     ], style = {'display': 'flex', 'flex-direction': 'row', 'flex-wrap': 'wrap', 'justify-content': 'left', 'width': '100%', 'height': '5vh', 'align-items': 'center', 'margin-top': '1vh'})
    
    # svg2 = visualizer.highlightScores(mol1, scores)




    # if (smiles2 is not None and len(smiles2) > 0):
    #     return dash_svg(svg1), dash_svg(svg2), "Success", base64.b64encode(pickle.dumps(siteLocator)).decode(), stats 
    # return dash_svg(svg1), dash_svg(svg2), "Success", base64.b64encode(pickle.dumps(siteLocator)).decode(), stats