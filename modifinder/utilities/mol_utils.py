"""
GNPS Utils - Molecule Utils
---------------------------
This file contains utility functions around molecules

Author: Shahneh
"""

import copy
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdFMCS
from copy import deepcopy
from modifinder.utilities.network import *

def get_modification_nodes(mol1, mol2, in_mol1 = True):
    """
        Calculates the modification sites between two molecules when one molecule is a substructure of the other molecule
        Input:
            mol1: first molecule
            mol2: second molecule
            in_mol1: bool, if True, the modification sites are given in the mol1, if False, the modification sites are given in the mol2
        Output:
            list of modification sites
    """

    copy_mol1, copy_mol2 = _get_molecules(mol1, mol2)
    
    if copy_mol1.HasSubstructMatch(copy_mol2):
        return _calculateModificationSites(copy_mol1, copy_mol2, in_mol1)
    elif copy_mol2.HasSubstructMatch(copy_mol1):
        return _calculateModificationSites(copy_mol2, copy_mol1, not in_mol1)
    else:
        raise ValueError("One molecule should be a substructure of the other molecule")


def get_modification_edges(mol1, mol2, only_outward_edges = False):
    """
        Calculates the modification edges between two molecules when one molecule is a substructure of the other molecule
        Input:
            mol1: first molecule
            mol2: second molecule
            only_outward_edges: bool, if True, only the modification edges that go from atoms in the substructure to atoms outside the substructure are returned
        Output:
            list of the indices of the modification edges in the parent molecule
    """

    copy_mol1, copy_mol2 = _get_molecules(mol1, mol2)
    
    if not mol1.HasSubstructMatch(mol2):
        if mol2.HasSubstructMatch(mol1):
            copy_mol1, copy_mol2 = copy_mol2, copy_mol1
        else:
            raise ValueError("One molecule should be a substructure of the other molecule")
    
    matches = _find_minimal_modification_edges_match(copy_mol1, copy_mol2)
    modificationEdgesOutward, modificationEdgesInside = _get_edge_modifications(copy_mol1, copy_mol2, matches)
    if only_outward_edges:
        return modificationEdgesOutward
    else:
        return modificationEdgesOutward + modificationEdgesInside


def get_edit_distance(mol1, mol2):
    """
        Calculates the edit distance between mol1 and mol2.
        Input:
            mol1: first molecule
            mol2: second molecule
        Output:
            edit_distance: edit distance between mol1 and mol2
    """
    copy_mol1, copy_mol2 = _get_molecules(mol1, mol2)
    mcs1 = rdFMCS.FindMCS([copy_mol1, copy_mol2])
    mcs_mol = Chem.MolFromSmarts(mcs1.smartsString)
    dist1 = get_modification_edges(copy_mol1, mcs_mol)
    dist2 = get_modification_edges(copy_mol2, mcs_mol)
    return len(dist1) + len(dist2)


def get_transition(input1, input2):
    """
        Calculates the transition between mol1 and mol2.
        Input:
            input1: first molecule
            input2: second molecule
        Output:
            result: a dictionary with the following keys:
                'merged_mol': the merged molecule
                'common_bonds': the common bonds between mol1 and mol2
                'common_atoms': the common atoms between mol1 and mol2
                'removed_atoms': the removed atoms from mol1
                'added_atoms': the added atoms from mol2
                'added_edges_inside': the added edges inside the common substructure
                'added_edges_bridge': the added edges between the common substructure and the added atoms
                'removed_edges_inside': the removed edges inside the common substructure
                'removed_edges_bridge': the removed edges between the common substructure and the removed
    """
    mol1, mol2 = _get_molecules(input1, input2)
    copy_mol1, copy_mol2 = Chem.Mol(mol1), Chem.Mol(mol2)
    # finding the maximum common substructure
    mcs1 = rdFMCS.FindMCS([mol1, mol2])
    mcs_mol = Chem.MolFromSmarts(mcs1.smartsString)

    # finding the mapping between the atoms of the common substructure
    mol1_indices = _find_minimal_modification_edges_match(mol1, mcs_mol)
    mol2_indices = _find_minimal_modification_edges_match(mol2, mcs_mol)

    map_mcs_from_mol1_to_mol2 = dict()
    map_mcs_from_mol2_to_mol1 = dict()
    for i in range(mcs_mol.GetNumAtoms()):
        map_mcs_from_mol1_to_mol2[mol1_indices[i]] = mol2_indices[i]
        map_mcs_from_mol2_to_mol1[mol2_indices[i]] = mol1_indices[i]

    # calculating the removed atoms from mol1
    removed_atoms = list(set(range(mol1.GetNumAtoms())) - set(mol1_indices))
    removed_edges_bridge, removed_edges_inside = _get_edge_modifications(mol1, mcs_mol, mol1_indices)

    # calculating the added atoms from mol2
    added_atoms = set(range(mol2.GetNumAtoms())) - set(mol2_indices)
    added_edges_bridge, added_edges_inside = _get_edge_modifications(mol2, mcs_mol, mol2_indices)

    # creating the merged molecule
    editable_mol = Chem.EditableMol(copy_mol1)
    ## adding the atoms from mol2
    n = mol1.GetNumAtoms()
    map_old_mol2_to_added = dict() # mapping from the old mol2 indices to the new indices in the merged molecule (for the added atoms)
    for atom in added_atoms:
        editable_mol.AddAtom(mol2.GetAtomWithIdx(atom))
        map_old_mol2_to_added[atom] = n
        n += 1
    for atom in mol2_indices:
        map_old_mol2_to_added[atom] = map_mcs_from_mol2_to_mol1[atom]
    print("map_old_mol2_to_added", map_old_mol2_to_added)
    ## adding the bonds from mol2
    updated_added_edges_inside = []
    for bond in added_edges_inside:
        print(bond)
        editable_mol.AddBond(map_old_mol2_to_added[bond[0]], map_old_mol2_to_added[bond[1]], order=copy_mol2.GetBondBetweenAtoms(bond[0], bond[1]).GetBondType())
        updated_added_edges_inside.append((map_old_mol2_to_added[bond[0]], map_old_mol2_to_added[bond[1]]))
    updated_added_edges_bridge = []
    for bond in added_edges_bridge:
        if bond[0] in added_atoms:
            index1 = map_old_mol2_to_added[bond[0]]
            index2 = map_mcs_from_mol2_to_mol1[bond[1]]
        else:
            index1 = map_mcs_from_mol2_to_mol1[bond[0]]
            index2 = map_old_mol2_to_added[bond[1]]
        updated_added_edges_bridge.append((index1, index2))
        editable_mol.AddBond(index1, index2, order=copy_mol2.GetBondBetweenAtoms(bond[0], bond[1]).GetBondType())
    ## adding the common bonds
    common_bonds = []
    for bond in mol1.GetBonds():
        pair = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        if pair in removed_edges_inside or pair in removed_edges_bridge:
            continue
        pair2 = (bond.GetEndAtomIdx(), bond.GetBeginAtomIdx())
        if pair2 in removed_edges_inside or pair2 in removed_edges_bridge:
            continue
        common_bonds.append(pair)

    
    common_atoms = mol1_indices
    mol1_atoms = removed_atoms
    mol2_atoms = [map_old_mol2_to_added[i] for i in added_atoms]

    # create a structure for result
    result = {
        'merged_mol': editable_mol.GetMol(),
        'common_bonds': common_bonds,
        'common_atoms': common_atoms,
        'removed_atoms': mol1_atoms,
        'added_atoms': mol2_atoms,
        'added_edges_inside': updated_added_edges_inside,
        'added_edges_bridge': updated_added_edges_bridge,
        'removed_edges_inside': removed_edges_inside,
        'removed_edges_bridge': removed_edges_bridge
    }

    return result


def get_modification_graph(main_struct, sub_struct):
    """
        Calculates the substructure difference between main_struct and sub_struct when there is exactly one modification edge,
          if there are multiple modification edges, one of the modifications is returned randomly.
        Input:
            main_struct: main molecule
            sub_struct: substructure molecule
        Output:
            frag: modified fragment mol
            index_in_frag: index of the modification atom in the fragment
            bondType: bond type of the modification bond
    """
    main_struct, sub_struct = _get_molecules(main_struct, sub_struct)
    if not main_struct.HasSubstructMatch(sub_struct):
        if sub_struct.HasSubstructMatch(main_struct):
            return get_modification_graph(sub_struct, main_struct)
        else:
            raise ValueError("One molecule should be a substructure of the other molecule")

    atoms_of_substructure = _find_minimal_modification_edges_match(main_struct, sub_struct)
    modificationEdgesOutward, modificationEdgesInside = _get_edge_modifications(main_struct, sub_struct, atoms_of_substructure)
    for bond in modificationEdgesOutward:
        frag_modif_atom = bond.GetBeginAtomIdx()
        bondType = bond.GetBondType()
        if bond.GetBeginAtomIdx() in atoms_of_substructure:
            frag_modif_atom = bond.GetEndAtomIdx()
        break

    # create a copy of the main structure
    main_struct_copy = Chem.Mol(main_struct)
    main_struct_copy.GetAtomWithIdx(frag_modif_atom).SetProp("atomNote", "modification")
    emol = Chem.EditableMol(main_struct_copy)

    for atom in reversed(range(main_struct_copy.GetNumAtoms())):
        if atom in atoms_of_substructure:
            emol.RemoveAtom(atom)
    
    frag = emol.GetMol()
    # get the atom index of the modification atom in the fragment
    for i in frag.GetAtoms():
        if i.HasProp("atomNote"):
            index_in_frag = i.GetIdx()
            i.ClearProp("atomNote")
            break
    
    return frag, index_in_frag, bondType


def attach_mols(main_mol, attachment_mol, attach_location_main, attach_location_attachment, bond_type):
    """
        Attaches the attachment structure to main molecule at the attach_location_main and attach_location_attachment with bond_type.
        Input:
            main_mol: rdkit molecule of the main molecule
            attachment_mol: rdkit molecule of the attachment molecule
            attach_location_main: the index of the atom in the main molecule where the attachment should be done
            attach_location_attachment: the index of the atom in the attachment molecule where the attachment should be done
            bond_type: the type of the bond between the main molecule and the attachment molecule
        Output:
            new_mol: the new molecule after attachment
    """

    new_mol = Chem.CombineMols(main_mol, attachment_mol)
    new_mol = Chem.EditableMol(new_mol)
    new_mol.AddBond(attach_location_main, attach_location_attachment + main_mol.GetNumAtoms(), bond_type)
    new_mol = new_mol.GetMol()
    return new_mol


def generate_possible_stuctures(main_struct, sub_struct):
    """
        Generates all possible structures after attaching sub_struct to main_struct.
        Input:
            main_struct: main molecule
            sub_struct: substructure molecule
        Output:
            list of possible_structures: all possible structures after attachment with the index of the atom
    """
    main_struct, sub_struct = _get_molecules(main_struct, sub_struct)
    if not main_struct.HasSubstructMatch(sub_struct):
        if sub_struct.HasSubstructMatch(main_struct):
            return generate_possible_stuctures(sub_struct, main_struct)
        else:
            raise ValueError("One molecule should be a substructure of the other molecule")

    frag, index_in_frag, bondType = get_modification_graph(main_struct, sub_struct)

    structs = []
    for atom in sub_struct.GetAtoms():
        temp_struct = _attach_struct_try(sub_struct, frag, atom.GetIdx(), index_in_frag, bondType)
        if temp_struct is None:
            continue
        structs.append((atom.GetIdx(), temp_struct))
    
    return structs
    

def _attach_struct_try(main_struct, frag, main_location, frag_location, bond_type):
    mol = attach_mols(main_struct, frag, main_location, frag_location, bond_type)
    try:
        Chem.SanitizeMol(mol)
        return mol
    except:
        for bond in [Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]:
            if bond != bond_type:
                try:
                    mol = attach_mols(main_struct, frag, main_location, frag_location, bond)
                    Chem.SanitizeMol(mol)
                    return mol
                except:
                    continue
    return None


def _find_minimal_modification_sites_match(mol, substructure):
    """
        Finds the matching that results in the minimum number of modification sites, if multiple matches are found.
        If multiple results are found, one random result is returned.
        Input:
            mol: first molecule
            substructure: substructure molecule
        Output:
            match array
    """
    matches = mol.GetSubstructMatches(substructure)
    if len(matches) == 1:
        return matches[0]
    else:
        min_modifications = float('inf')
        min_match = None
        for match in matches:
            modif_sites = set()
            for bond in mol.GetBonds():
                if bond.GetBeginAtomIdx() in match and ~(bond.GetEndAtomIdx() in match):
                    modif_sites.add(bond.GetBeginAtomIdx())
                if bond.GetEndAtomIdx() in match and ~(bond.GetBeginAtomIdx() in match):
                    modif_sites.add(bond.GetEndAtomIdx())
            if len(modif_sites) < min_modifications:
                min_modifications = len(modif_sites)
                min_match = match
        return min_match


def _find_minimal_modification_edges_match(mol, substructure):
    """
        Finds the matching that results in the minimum number of modification edges, if multiple matches are found.
        If multiple results are found, one random result is returned.
        Input:
            mol: first molecule
            substructure: substructure molecule
        Output:
            match array
    """
    if not mol.HasSubstructMatch(substructure):
        raise ValueError("The substructure is not a substructure of the molecule")
    matches = mol.GetSubstructMatches(substructure)
    if len(matches) == 1:
        return matches[0]
    else:
        min_modifications = float('inf')
        min_match = None
        for match in matches:
            modificationEdgesOutward, modificationEdgesInside = _get_edge_modifications(mol, substructure, match)
            num_modifications = len(modificationEdgesOutward) + len(modificationEdgesInside)
            if num_modifications < min_modifications:
                min_modifications = num_modifications
                min_match = match
        return min_match


def _get_edge_modifications(mol, substructure, match):
    """
        Calculates the modification edges to get mol from substructure given the mapping of the atoms
        Input:
            mol: first molecule
            substructure: substructure molecule
            match: the match array, mapping atoms from substructure to mol
        Output:
            count: modification edges
    """
    reverse_match = {v: k for k, v in enumerate(match)}
    intersect = set(match)
    modificationEdgesOutward = []
    modificationEdgesInside = []
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in intersect or bond.GetEndAtomIdx() in intersect:
            if bond.GetBeginAtomIdx() in intersect and bond.GetEndAtomIdx() in intersect:
                a1 = reverse_match[bond.GetBeginAtomIdx()]
                a2 = reverse_match[bond.GetEndAtomIdx()]
                if substructure.GetBondBetweenAtoms(a1, a2) is None:
                    modificationEdgesInside.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
            else:
                modificationEdgesOutward.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        else:
            continue
    
    return modificationEdgesOutward, modificationEdgesInside


def _calculateModificationSites(mol, substructure, inParent = True):
    """
        Calculates the number of modification sites to get mol from substructure (works when one moleculte is a substructure of the other molecule)
        if multiple matches are found, the modification sites that result in the minimum number of modifications are returned
        Input:
            mol1: first molecule
            substructure: substructure molecule
            inParent: bool, if True, the modification sites are given in the parent molecule, if False, the modification sites are given in the substructure
        Output:
            count: modification sites
    """

    matches = _find_minimal_modification_edges_match(mol, substructure)
    intersect = set(matches)

    
    modificationSites = []
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in intersect:
            for tempAtom in intersect:
                if mol.GetBondBetweenAtoms(atom.GetIdx(), tempAtom) is not None:
                    modificationSites.append(tempAtom)

    if inParent:
        return modificationSites
    else:
        res = []
        if type(matches[0]) is tuple:
            # multiple matches (substructure is happening multiple times in the parent molecule)
            for match in matches:
                subMatches = list(match)
                temp_res = []
                for atom in modificationSites:
                    idx = subMatches.index(atom) # based on rdkit library, the value in matches array are sorted by index
                    temp_res.append(idx)
                res.append(temp_res)
        else:
            for atom in modificationSites:
                subMatches = list(matches)
                idx = subMatches.index(atom) # based on rdkit library, the value in matches array are sorted by index
                res.append(idx)
        return res


def _get_molecules(mol1, mol2):
    copy_mol1 = copy.deepcopy(mol1)
    copy_mol2 = copy.deepcopy(mol2)
    copy_mol1 = _get_molecule(copy_mol1)
    copy_mol2 = _get_molecule(copy_mol2)
    if copy_mol1 is None or copy_mol2 is None:
        raise ValueError("The molecules are not valid")
    return copy_mol1, copy_mol2


def _get_molecule(identifier) -> Chem.Mol:
    """
    Get the rdkit molecule from smiles, inchi or GNPS identifier
    :param identifier: rdkit Mol or str - GNPS identifier (USI or Accession) or SMILES or InChI or MolBlock or SMARTS
    :return: rdkit molecule or None
    """
    if isinstance(identifier, Chem.Mol):
        return identifier
    try:
        molecule = Chem.MolFromSmiles(identifier)
        if molecule:
            return molecule
    except:
        pass

    try:
        molecule = Chem.MolFromInchi(identifier)
        if molecule:
            return molecule
    except:
        pass

    try:
        molecule = Chem.MolFromMolBlock(identifier)
        if molecule:
            return molecule
    except:
        pass

    try:
        molecule = Chem.MolFromSmarts(identifier)
        if molecule:
            return molecule
    except:
        pass

    try:
        data = get_data(identifier)
        smiles = data.get("Smiles", None)
        if smiles:
            molecule = Chem.MolFromSmiles(smiles)
            if molecule:
                return molecule
    except:
        raise ValueError("The identifier is not valid or it is not supported")
