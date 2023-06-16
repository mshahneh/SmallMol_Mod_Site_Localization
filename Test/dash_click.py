from dash import Dash, html, dcc, Input, Output, State
import rdkit.Chem as Chem
import rdkit.Chem.Draw as Draw
import dash_bootstrap_components as dbc

app = Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

if __name__ == '__main__':
    mol = Chem.MolFromSmiles("C1=CC=C(C=C1)C(=O)O")
    def mol_with_atom_index(mol):
        for atom in mol.GetAtoms():
            atom.SetProp("atomLabel", atom.GetSymbol())
        return mol
    mol = mol_with_atom_index(mol)
    d2d = Chem.Draw.MolDraw2DSVG(300,400)
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    molSVG = d2d.GetDrawingText()
    # remove the xml header
    molSVG = molSVG[molSVG.find("<svg"):]
    
    # add viewBox to svg element
    molSVG = molSVG.replace("viewBox='0 0 300 400'", "viewBox='0 0 300 300'")

    
    molSVG = molSVG.replace("width='300px'", "width='100%'")
    molSVG = molSVG.replace("height='400px'", "height='auto'")

    print(molSVG[0:350])
#     molSVG += '''<script>
#     var iframe = document.currentScript.previousElementSibling;
#     var svg = iframe.contentDocument.getElementById('svg');
    
#     function updateViewBox() {
#         var width = iframe.clientWidth;
#         var height = iframe.clientHeight;
#         svg.setAttribute('viewBox', '0 0 ' + width + ' ' + height);
#     }
    
#     window.addEventListener('resize', updateViewBox);
#     updateViewBox();
# </script>'''

    svg_data = '''
<svg xmlns="http://www.w3.org/2000/svg" width="200" height="200">
    <circle cx="100" cy="100" r="80" fill="blue" />
</svg>
'''

    # show svg in iframe, force it to scale into the container
    app.layout = html.Div([
        html.Iframe(
            id='mol',
            srcDoc=molSVG,
            style={
                'width': '10vw',
                'height': '40vh',
                'border': 'none',
            }
        ),
        html.Div(id='clicked_atom', style={'margin-top': '1em'}),
    ])

    # if a path in mol is clicked, show it in the div
    @app.callback(
        Output('clicked_atom', 'children'),
        [Input('mol', 'n_clicks')],
        [State('mol', 'clickData')])
    def display_click_data(n_clicks, clickData):
        if clickData is None:
            return 'No atom clicked.'
        atom = clickData['points'][0]['text']
        return f'Atom {atom} was clicked.'
    

    


    app.run_server(debug = True, port=8051)