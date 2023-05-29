import xlsxwriter
import os
import utils as utils
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
from io import BytesIO
import cairosvg
import pandas as pd
import urllib.parse
import matplotlib.pyplot as plt

def make_url(base_url = "http://localhost:8050/", USI1=None, USI2=None, SMILES1=None, SMILES2=None, args=None):
    query_params = {k: v for k, v in {
        "USI1": USI1,
        "USI2": USI2,
        "SMILES1": SMILES1,
        "SMILES2": SMILES2
    }.items() if v is not None}
    
    if args is not None:
        query_params.update(args)

    url = base_url + "?" + urllib.parse.urlencode(query_params)
    return url

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
            d2d = Chem.Draw.MolDraw2DSVG(1250,1200)
            d2d.DrawMolecule(mol)
        else:
            hitAtoms, hitBonds = utils.getHitAtomsAndBonds(mol, substructure)
            if highlightModificationSites:
                modifications = utils.calculateModificationSites(mol, substructure)
                d2d = Chem.Draw.MolDraw2DSVG(250, 200)
                colors = dict()
                for hitAtom in hitAtoms[0]:
                    if (hitAtom in modifications):
                        colors[hitAtom] = (0.3, 0.6, 1)
                    else:
                        colors[hitAtom] = (1,0.5,0.5)
                d2d.DrawMolecule(mol,highlightAtoms=hitAtoms[0], highlightBonds=hitBonds[0], highlightAtomColors=colors)
            else:
                d2d = Chem.Draw.MolDraw2DSVG(1250,1200)
                d2d.DrawMolecule(mol,highlightAtoms=hitAtoms[0], highlightBonds=hitBonds[0])
    else:
        d2d = Chem.Draw.MolDraw2DSVG(1250,1200)
        d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    return d2d.GetDrawingText()

def highlightScores(mol, scores):
    vals = [x/max(scores) for x in scores]
    d2d = Draw.MolDraw2DSVG(1250,1200)
    colors = dict()
    for i in range(0, mol.GetNumAtoms()):
        colors[i] = (1-(vals[i]**2), 1-(vals[i]**2), 1-(vals[i]**2))
    d2d.DrawMolecule(mol, highlightAtoms=list(range(mol.GetNumAtoms())), highlightAtomColors=colors, highlightBonds=[])
    d2d.FinishDrawing()
    return d2d.GetDrawingText()

def highlightMolIndices(mol, hitAtoms):
    d2d = Draw.MolDraw2DSVG(1250,1200)
    hitBonds = []
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in hitAtoms and bond.GetEndAtomIdx() in hitAtoms:
            hitBonds.append(bond.GetIdx())
    d2d = Chem.Draw.MolDraw2DSVG(1250,1200)
    d2d.DrawMolecule(mol,highlightAtoms=hitAtoms, highlightBonds=hitBonds)
    d2d.FinishDrawing()
    return d2d.GetDrawingText()

def draw_alignment(peaks1, peaks2, matched_peaks, shift = 0.1, show_text = False, show_lines = True, scale = 1, ax = None):

    if ax is None:
        fig, ax = plt.subplots(figsize=(20*scale, 10*scale))

    shifted, unshifted = utils.separateShifted(matched_peaks, peaks1, peaks2)
    max_y = max([_[1] for _ in peaks1])
    peaks1 = [(peaks1[i][0], peaks1[i][1]/max_y) for i in range(len(peaks1))]

    max_y = max([_[1] for _ in peaks2])
    peaks2 = [(peaks2[i][0], peaks2[i][1]/max_y) for i in range(len(peaks2))]

    max_x = max([_[0] for _ in peaks1])
    max_x = max([max_x, max([_[0] for _ in peaks2])])

    ax.set_xlim(0, max_x+5)

    print('debugging', peaks2)
    
    molShifted = [i[0] for i in shifted]
    molUnshifted = [i[0] for i in unshifted]
    modifiedShifted = [i[1] for i in shifted]
    modifiedUnshifted = [i[1] for i in unshifted]

    # draw peaks in site locator main compound
    # make shifted peaks red
    # make unshifted peaks blue

    x = [_[0] for _ in peaks1 ]
    y = [ _[1] for _ in peaks1 ]
    shift_array = [shift for _ in peaks1 ]
    ax.bar(x, shift_array, width=0.5*scale, color="white", alpha = 0)
    ax.bar(x, y, width=0.5*scale, color="gray", bottom=shift_array)

    x_shifted = [peaks1[_][0] for _ in molShifted]
    y_shifted = [ peaks1[_][1] for _ in molShifted]
    shift_array = [shift for _ in molShifted]
    ax.bar(x_shifted, shift_array, width=0.5*scale, color="white", alpha = 0)
    ax.bar(x_shifted, y_shifted, width=0.5*scale, color="red", bottom=shift_array)

    x_unshifted = [peaks1[_][0] for _ in molUnshifted]
    y_unshifted = [ peaks1[_][1] for _ in molUnshifted]
    shift_array = [shift for _ in molUnshifted]
    ax.bar(x_unshifted, shift_array, width=0.5*scale, color="white", alpha = 0)
    ax.bar(x_unshifted, y_unshifted, width=0.5*scale, color="blue", bottom=shift_array)

    #plot modified peaks as reversed
    x = [_[0] for _ in peaks2]
    y = [-_[1] for _ in peaks2]
    ax.bar(x, y, width=0.5*scale, color="gray")

    x_shifted = [peaks2[_][0] for _ in modifiedShifted]
    y_shifted = [-peaks2[_][1] for _ in modifiedShifted]
    ax.bar(x_shifted, y_shifted, width=0.5*scale, color="red")

    x_unshifted = [peaks2[_][0] for _ in modifiedUnshifted]
    y_unshifted = [-peaks2[_][1] for _ in modifiedUnshifted]
    ax.bar(x_unshifted, y_unshifted, width=0.5*scale, color="blue")

    if show_lines:
        for peak in shifted:
            x1 = peaks1[peak[0]][0]
            x2 = peaks2[peak[1]][0]
            # draw line with text in the middle
            ax.plot([x1, x2], [shift, 0], color="red", linewidth=0.5*scale, linestyle="--")
            val = abs(peaks1[peak[0]][0] - peaks2[peak[1]][0])
            val = round(val, 2)
            if show_text:
                ax.text((x1 + x2) / 2, shift/2, str(val), fontsize=10, horizontalalignment='center')

        for peak in unshifted:
            x1 = peaks1[peak[0]][0]
            x2 = peaks2[peak[1]][0]
            ax.plot([x1, x2], [shift, 0], color="blue", linewidth=0.5*scale, linestyle="--")

    #draw horizontal line
    ax.plot([5, max_x], [shift, shift], color="gray", linewidth=0.5*scale, linestyle="-")
    ax.plot([5, max_x], [0, 0], color="gray", linewidth=0.5*scale, linestyle="-")

    # custom y axis ticks
    y_ticks1 = [i/10 + shift for i in range(0, 11, 2)]
    y_ticks2 = [-i/10 for i in range(0, 11, 2)]
    # reverse y ticks2
    y_ticks2 = y_ticks2[::-1]
    y_ticks =  y_ticks2 + y_ticks1 
    print(y_ticks)
    y_ticks_labels1 = [i for i in range(0, 110, 20)]
    y_ticks_labels2 = [i for i in range(0, 110, 20)]
    # reverse y ticks2
    y_ticks_labels2 = y_ticks_labels2[::-1]
    y_ticks_labels = y_ticks_labels2 + y_ticks_labels1
    print(y_ticks_labels)
    ax.set_yticks(y_ticks, y_ticks_labels)

    # set font size
    ax.tick_params(axis='both', which='major', labelsize=20*scale)

    # add legend
    ax.plot([], [], color="red", linewidth=2*scale, linestyle="-", label="Shifted Matched Peaks")
    ax.plot([], [], color="blue", linewidth=2*scale, linestyle="-", label="Unshifted Matched Peaks")
    ax.plot([], [], color="gray", linewidth=2*scale, linestyle="-", label="Unmatched Peaks")
    ax.legend(loc='upper left')
    # legend font size
    ax.legend(prop={'size': 20*scale})

    return ax