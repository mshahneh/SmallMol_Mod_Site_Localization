import xlsxwriter
import os
import utils as utils
from rdkit import Chem
from PIL import Image
from io import BytesIO
import cairosvg

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

def molToSVG(mol, substructure=None, highlightModificationSites=False):
    """
    Converts a molecule to an SVG image.
    """
    if substructure is not None:
        hitAtoms, hitBonds = utils.getHitAtomsAndBonds(mol, substructure)
        if highlightModificationSites:
            modifications = list(utils.calculateModificationSites(mol, substructure))
            d2d = Chem.Draw.MolDraw2DSVG(250,200)
            colors = dict()
            for hitAtom in hitAtoms[0]:
                if (hitAtom in modifications):
                    colors[hitAtom] = (0.5,0.5,1)
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