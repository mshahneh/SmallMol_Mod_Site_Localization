from rdkit import Chem
from rdkit.Chem import Draw
import traceback

def getHitAtomsAndBonds(mol, substructure):
    """
        Returns the atoms and bonds that are hit by the substructure.
        Input:
            mol: molecule
            substructure: substructure molecule
        Output:
            hitAtoms: list of hit atoms
            hitBonds: list of hit bonds
    """
    
    matches = mol.GetSubstructMatches(substructure)
    hitAtoms = []
    hitBonds = []
    if len(matches) == 0:
        return hitAtoms, hitBonds
    
    if type(matches[0]) is tuple:
        for match in matches:
            tempHitAtoms = set()
            tempHitBonds = set()
            for atom in match:
                tempHitAtoms.add(atom)
            for bond in mol.GetBonds():
                if bond.GetBeginAtomIdx() in tempHitAtoms and bond.GetEndAtomIdx() in tempHitAtoms:
                    tempHitBonds.add(bond.GetIdx())
            if len(tempHitAtoms) > 0:
                hitAtoms.append(list(tempHitAtoms))
                hitBonds.append(list(tempHitBonds))
    else:
        hitAtoms.append([])
        hitBonds.append([])
        for atom in matches:
            hitAtoms[0].append(atom)
        for bond in mol.GetBonds():
            if bond.GetBeginAtomIdx() in hitAtoms[0] and bond.GetEndAtomIdx() in hitAtoms[0]:
                hitBonds[0].append(bond.GetIdx())

    return hitAtoms, hitBonds

def molToSVG(mol, substructure=None, highlightModificationSites=False):
    """
    Converts a molecule to an SVG image.
    """
    if substructure is not None:
        if (type(substructure) == str):
            substructure = Chem.MolFromSmarts(substructure)
        if not mol.HasSubstructMatch(substructure):
            d2d = Draw.MolDraw2DSVG(250,200)
            d2d.DrawMolecule(mol)
        else:
            hitAtoms, hitBonds = getHitAtomsAndBonds(mol, substructure)
            if highlightModificationSites:
                modifications = list(calculateModificationSites(mol, substructure))
                d2d = Draw.MolDraw2DSVG(250, 200)
                colors = dict()
                for hitAtom in hitAtoms[0]:
                    if (hitAtom in modifications):
                        colors[hitAtom] = (0.3, 0.6, 1)
                    else:
                        colors[hitAtom] = (1,0.5,0.5)
                d2d.DrawMolecule(mol,highlightAtoms=hitAtoms[0], highlightBonds=hitBonds[0], highlightAtomColors=colors)
            else:
                d2d = Draw.MolDraw2DSVG(250,200)
                d2d.DrawMolecule(mol,highlightAtoms=hitAtoms[0], highlightBonds=hitBonds[0])
    else:
        d2d = Draw.MolDraw2DSVG(250,200)
        d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    return d2d.GetDrawingText()

def calculateModificationSites(mol, substructure, inParent = True):
    """
        Calculates the number of modification sites to get mol from substructure.
        Input:
            mol1: first molecule
            substructure: substructure molecule
        Output:
            count: modification sites
    """

    matches = mol.GetSubstructMatch(substructure)
    intersect = set(matches)

    
    modificationSites = set()
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in intersect:
            for tempAtom in intersect:
                if mol.GetBondBetweenAtoms(atom.GetIdx(), tempAtom) is not None:
                    modificationSites.add(tempAtom)

    if inParent:
        return modificationSites
    else:
        res = set()
        for atom in modificationSites:
            if type(matches[0]) is tuple:
                for match in matches:
                    subMatches = list(match)
                    idx = subMatches.index(atom)
                    res.add(idx)
            else:
                subMatches = list(matches)
                idx = subMatches.index(atom)
                res.add(idx)
        return res

def locate_modifications(sp1, sp2, output, output_format="svg"):
    output_path = output + "." + output_format
    try:
        m1 = Chem.MolFromSmiles(sp1)
        m2 = Chem.MolFromSmiles(sp2)
        if m1.HasSubstructMatch(m2):
            svg = molToSVG(m1, m2, True)
        elif m2.HasSubstructMatch(m1):
            svg = molToSVG(m2, m1, True)
        else:
            print("Does not qualify for modification site calculation")
            return
        
        if output_format == "png":
            import cairosvg
            cairosvg.svg2png(bytestring=svg, write_to=output_path)
        else:    
            with open(output_path, "w") as f:
                f.write(svg)
    except:
        print("Exception occurred")
        traceback.print_exc()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='library downloader')
    parser.add_argument('--sp1', type=str, default="C1=C(C=C(C=C1N)O)C(=O)O", help='Spectra 1')
    parser.add_argument('--sp2', type=str, default="OCC1=CC=CC=C1", help='Spectra 2')
    parser.add_argument('--output', type=str, default="output", help='output file name (without extension))')
    parser.add_argument('--output_format', type=str, default="svg", help='output format (svg, png)')
    args = parser.parse_args()
    locate_modifications(**vars(args))