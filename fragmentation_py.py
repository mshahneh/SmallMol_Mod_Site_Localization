# 2014.10.20 12:28:02 CEST
import numpy
import pars as pars
import os
from rdkit import Chem

class FragmentEngine(object):

    def __init__(self, mol, max_broken_bonds, max_water_losses, ionisation_mode, skip_fragmentation, molcharge):
        try:
            self.mol = Chem.MolFromMolBlock(str(mol))
            self.accept = True
            self.natoms = len(self.mol.GetAtoms())
        except:
            self.accept = False
            return 
        self.max_broken_bonds = max_broken_bonds
        self.max_water_losses = max_water_losses
        self.ionisation_mode = ionisation_mode
        self.skip_fragmentation = skip_fragmentation
        self.molcharge = molcharge
        self.atom_masses = []
        self.neutral_loss_atoms = []
        self.bonded_atoms = []
        self.bonds = set([])
        self.bondscore = {}
        self.new_fragment = 0
        self.template_fragment = 0
        self.fragment_masses = ((max_broken_bonds + max_water_losses) * 2 + 1) * [0]
        self.fragment_info = [[0, 0, 0]]
        self.avg_score = None
        for x in self.mol.GetAtoms():
            self.bonded_atoms.append([])
            self.atom_masses.append(x.GetMass())
            if x.GetSymbol() == 'O' and x.GetTotalNumHs() == 1 and x.GetTotalDegree() == 1:
                self.neutral_loss_atoms.append(x.GetIdx())
            if x.GetSymbol() == 'N' and x.GetTotalNumHs() == 2 and x.GetTotalDegree() == 1:
                self.neutral_loss_atoms.append(x.GetIdx())

        for x in range(self.mol.GetNumBonds()):
            (a1, a2,) = self.mol.GetBondWithIdx(x).GetBeginAtomIdx(), self.mol.GetBondWithIdx(x).GetEndAtomIdx()
            self.bonded_atoms[a1].append(a2)
            self.bonded_atoms[a2].append(a1)
            bond = 1 << a1 | 1 << a2
            bondscore = pars.typew[self.mol.GetBondWithIdx(x).GetBondType().name] * pars.heterow[(self.mol.GetAtomWithIdx(a1).GetSymbol() != 'C' or self.mol.GetAtomWithIdx(a2).GetSymbol() != 'C')]
            self.bonds.add(bond)
            self.bondscore[bond] = bondscore




    def extend(self, atom):
        for a in self.bonded_atoms[atom]:
            atombit = 1 << a
            if atombit & self.template_fragment and not atombit & self.new_fragment:
                self.new_fragment = self.new_fragment | atombit
                self.extend(a)




    def generate_fragments(self):
        frag = (1 << self.natoms) - 1
        all_fragments = set([frag])
        total_fragments = set([frag])
        current_fragments = set([frag])
        new_fragments = set([frag])
        self.add_fragment(frag, self.calc_fragment_mass(frag), 0, 0)
        if self.skip_fragmentation:
            self.convert_fragments_table()
            return len(self.fragment_info)
        for step in range(self.max_broken_bonds):
            for fragment in current_fragments:
                for atom in range(self.natoms):
                    if 1 << atom & fragment:
                        self.template_fragment = fragment ^ 1 << atom
                        list_ext_atoms = set([])
                        extended_fragments = set([])
                        for a in self.bonded_atoms[atom]:
                            if 1 << a & self.template_fragment:
                                list_ext_atoms.add(a)

                        if len(list_ext_atoms) == 1:
                            extended_fragments.add(self.template_fragment)
                        else:
                            for a in list_ext_atoms:
                                for frag in extended_fragments:
                                    if 1 << a & frag:
                                        break
                                else:
                                    self.new_fragment = 1 << a
                                    self.extend(a)
                                    extended_fragments.add(self.new_fragment)


                        for frag in extended_fragments:
                            if frag not in all_fragments:
                                all_fragments.add(frag)
                                (bondbreaks, score,) = self.score_fragment(frag)
                                if bondbreaks <= self.max_broken_bonds and score < pars.missingfragmentpenalty + 5:
                                    new_fragments.add(frag)
                                    total_fragments.add(frag)
                                    self.add_fragment(frag, self.calc_fragment_mass(frag), score, bondbreaks)



            current_fragments = new_fragments
            new_fragments = set([])

        for step in range(self.max_water_losses):
            for fi in self.fragment_info:
                if fi[2] == self.max_broken_bonds + step:
                    fragment = fi[0]
                    for atom in self.neutral_loss_atoms:
                        if 1 << atom & fragment:
                            frag = fragment ^ 1 << atom
                            if frag not in total_fragments:
                                total_fragments.add(frag)
                                (bondbreaks, score,) = self.score_fragment(frag)
                                if score < pars.missingfragmentpenalty + 5:
                                    self.add_fragment(frag, self.calc_fragment_mass(frag), score, bondbreaks)



        self.convert_fragments_table()
        return len(self.fragment_info)



    def score_fragment(self, fragment):
        score = 0
        bondbreaks = 0
        for bond in self.bonds:
            if 0 < fragment & bond < bond:
                score += self.bondscore[bond]
                bondbreaks += 1

        if score == 0:
            print('score=0: ', end=' ')
            print(fragment, end=' ')
            print(bondbreaks)
        return (bondbreaks, score)



    def score_fragment_rel2parent(self, fragment, parent):
        score = 0
        for bond in self.bonds:
            if 0 < fragment & bond < bond & parent:
                score += self.bondscore[bond]

        return score



    def calc_fragment_mass(self, fragment):
        fragment_mass = 0.0
        for atom in range(self.natoms):
            if fragment & 1 << atom:
                fragment_mass += self.atom_masses[atom]

        return fragment_mass



    def add_fragment(self, fragment, fragmentmass, score, bondbreaks):
        mass_range = (self.max_broken_bonds + self.max_water_losses - bondbreaks) * [0] + list(numpy.arange(-bondbreaks + self.ionisation_mode * (1 - self.molcharge), bondbreaks + self.ionisation_mode * (1 - self.molcharge) + 1) * pars.Hmass + fragmentmass) + (self.max_broken_bonds + self.max_water_losses - bondbreaks) * [0]
        if bondbreaks == 0:
            mass_range[self.max_broken_bonds + self.max_water_losses - self.ionisation_mode] = fragmentmass
        self.fragment_masses += mass_range
        self.fragment_info.append([fragment, score, bondbreaks])



    def convert_fragments_table(self):
        self.fragment_masses_np = numpy.array(self.fragment_masses).reshape(len(self.fragment_info), (self.max_broken_bonds + self.max_water_losses) * 2 + 1)



    def calc_avg_score(self):
        self.avg_score = numpy.average(self.scores)



    def get_avg_score(self):
        return self.avg_score



    def find_fragments(self, mass, parent, precision, mz_precision_abs):
        result = numpy.where(numpy.where(self.fragment_masses_np < max(mass * precision, mass + mz_precision_abs), self.fragment_masses_np, 0) > min(mass / precision, mass - mz_precision_abs))
        fragment_set = []
        for i in range(len(result[0])):
            fid = result[0][i]
            fragment_set.append(self.fragment_info[fid] + [self.fragment_masses_np[fid][(self.max_broken_bonds + self.max_water_losses - self.ionisation_mode * (1 - self.molcharge))]] + [self.ionisation_mode * (1 - self.molcharge) + result[1][i] - self.max_broken_bonds - self.max_water_losses])

        return fragment_set



    def get_fragment_info(self, fragment, deltaH):
        atomstring = ''
        atomlist = []
        elements = dict([(e,0) for e in list(pars.mims.keys())])
        for atom in self.mol.GetAtoms():
            if 1 << atom.GetIdx() & fragment:
                atomstring += ',' + str(atom.GetIdx())
                atomlist.append(atom.GetIdx())
                elements[atom.GetSymbol()] += 1
                elements['H'] += atom.GetTotalNumHs()
        # for atom in range(self.natoms):
        #     if 1 << atom & fragment:
        #         atomstring += ',' + str(atom)
        #         atomlist.append(atom)
        #         elements[Chem.GetAtomSymbol(self.mol, atom)] += 1
        #         elements['H'] += Chem.GetAtomHs(self.mol, atom)

        formula = ''
        for el in list(pars.mims.keys()):
            nel = elements[el]
            if nel > 0:
                formula += el
            if nel > 1:
                formula += str(nel)

        #get inchi key of fragment
        fragmentSmiles = Chem.MolFragmentToSmiles(self.mol, atomlist)
        # fragmentMol = Chem.MolFromSmiles(fragmentSmiles)
        # fragmentInchi = Chem.MolToInchi(fragmentMol)

        return (atomstring,
         atomlist,
         formula,
         fragmentSmiles
        )



    def get_natoms(self):
        return self.natoms



    def accepted(self):
        return self.accept
