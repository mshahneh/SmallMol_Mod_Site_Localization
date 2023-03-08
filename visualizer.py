import xlsxwriter
import os
import utils as utils
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
from io import BytesIO
import cairosvg
import pandas as pd

def table_to_xlsx(data, path):
    """
    Converts a dictionary to an Excel file.
    """
    wb = xlsxwriter.Workbook(path)
    ws = wb.add_worksheet('Sheet 1')
    col = 0
    for key in data:
        if key != 'image':
            ws.write(0, col, key)
            ws.write(1, col, data[key])
            col += 1
        else:
            ws.write(0, col, key)
            # write svg to file
            with open("temp.svg", "w") as f:
                f.write(data[key])
            # svg image to bitmap
            cairosvg.svg2png(url="temp.svg", write_to="temp.png", parent_width=250, parent_height=200)
            # # convert png to bitmap
            # img = Image.open("temp.png")
            # image_parts = img.split()
            # r = image_parts[0]
            # g = image_parts[1]
            # b = image_parts[2]
            # img = Image.merge("RGB", (r, g, b))
            # fo = BytesIO()
            # img.save(fo, format='bmp')

            # ws.insert_bitmap_data(fo.getvalue(), col, 0)
            # img.close()
            ws.insert_image(1, col, "temp.png")


            col += 1
        
    wb.close()
    #delete temp files
    try:
        os.remove("temp.svg")
    except:
        pass
    try:
        os.remove("temp.png")
    except:
        pass

def table_to_tsv(data, path_to_tsv, path_to_img_folder):
    import uuid

    unique_key = str(uuid.uuid4())

    path_to_svg = os.path.join(path_to_img_folder, unique_key + ".svg")

    data["img_path"] = path_to_svg

    df = pd.DataFrame([data])
    df = df[['# matched peaks', '# unshifted peaks', '# shifted peaks', 'usi1', 'usi2', 'structure1', 'img_path']]
    
    # TODO: Adding structures and USIs

    df.to_csv(path_to_tsv, sep='\t', index=False)

    # Writing out the images to file
    with open(path_to_svg, "w") as f:
        f.write(data["image"])


def molToSVG(mol, substructure=None, highlightModificationSites=False):
    """
    Converts a molecule to an SVG image.
    """
    if substructure is not None:
        if (type(substructure) == str):
            substructure = Chem.MolFromSmarts(substructure)
        if not mol.HasSubstructMatch(substructure):
            d2d = Chem.Draw.MolDraw2DSVG(250,200)
            d2d.DrawMolecule(mol)
        else:
            hitAtoms, hitBonds = utils.getHitAtomsAndBonds(mol, substructure)
            if highlightModificationSites:
                modifications = list(utils.calculateModificationSites(mol, substructure))
                d2d = Chem.Draw.MolDraw2DSVG(250, 200)
                colors = dict()
                for hitAtom in hitAtoms[0]:
                    if (hitAtom in modifications):
                        colors[hitAtom] = (0.3, 0.6, 1)
                    else:
                        colors[hitAtom] = (1,0.5,0.5)
                d2d.DrawMolecule(mol,highlightAtoms=hitAtoms[0], highlightBonds=hitBonds[0], highlightAtomColors=colors)
            else:
                d2d = Chem.Draw.MolDraw2DSVG(250,200)
                d2d.DrawMolecule(mol,highlightAtoms=hitAtoms[0], highlightBonds=hitBonds[0])
    else:
        d2d = Chem.Draw.MolDraw2DSVG(250,200)
        d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    return d2d.GetDrawingText()

def highlightScores(mol, scores):
    d2d = Draw.MolDraw2DSVG(250,200)
    colors = dict()
    for i in range(0, mol.GetNumAtoms()):
        colors[i] = (1-(scores[i]**2), 1-(scores[i]**2), 1-(scores[i]**2))
    d2d.DrawMolecule(mol, highlightAtoms=list(range(mol.GetNumAtoms())), highlightAtomColors=colors)
    d2d.FinishDrawing()
    return d2d.GetDrawingText()