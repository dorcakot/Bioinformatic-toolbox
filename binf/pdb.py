from Bio.PDB import *
import math
import numpy as np
import matplotlib.pyplot as plt

POLAR = ['ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'HIS', 'LYS', 'SER', 'THR', 'TYR']

class PDBMolecule:
    """
    PDBMolecule represents a molecule stored in pdb format.
    """
    def __init__(self, file):
        """
        To initialize, use name of a file containing data in pdb format. This data is parsed and stored as 'structure'.
        :param file: name of a file containing pdb data
        """
        parser = PDBParser()
        self.structure = parser.get_structure(file.split('.')[0], file)
        self.model = self.structure.get_models()

    def get_chain(self, id):
        """
        Function to obtain an object representing a chain by its id.
        :param id: id of a chain
        :return: object representing a chain
        """
        for chain in self.structure.get_chains():
            if chain.id == id:
                return chain

    def get_residue(self, id):
        """
        Function to obtain an object representing a residue by its id.
        :param id: id of a residue
        :return: object representing a residue
        """
        for residue in self.structure.get_residues():
            if id in residue.id:
                return residue

    def get_atom(self, id, residue):
        """
        Function to obtain an object representing an atom specified by its residue and id.
        :param id: id of an atom
        :param residue: residue containing the atom
        :return: object representing an atom
        """
        for atom in residue:
            if id in atom.id:
                return atom

    def get_information(self):
        """
        Function to display information about number of models, chains, residues, atom within a pdb molecule.
        :return: information about a model
        """
        info = ""
        info += "number of models: "
        models = sum(1 for m in self.structure.get_models())
        info += str(models)
        info += "\nnumber of chains: "
        chains = sum(1 for m in self.structure.get_chains())
        info += str(chains)
        info += "\nnumber of residues: "
        residues = sum(1 for m in self.structure.get_residues())
        info += str(residues)
        info += "\nnumber of atoms: "
        atoms = sum(1 for m in self.structure.get_atoms())
        info += str(atoms)
        return info

    def get_width(self):
        """
        This function calculates width of the molecule by finding the greatest distance between two atoms.
        :return: widt of the molecule
        """
        width = 0
        atoms = self.structure.get_atoms()
        for atom in atoms:
            for a in atoms:
                a1 = atom.get_vector()
                a2 = a.get_vector()
                current = math.sqrt((a2[0] - a1[0])**2 + (a2[1] - a1[1])**2 + (a2[2] - a1[2])**2)
                if current > width:
                    width = current
        return round(width,3)

    def get_atoms_at_distance(self, hetatm, distance, round_digits=1):
        """
        Function serves to find all the atoms at given distance from a heteroatom.
        :param hetatm: object representing an atom
        :param distance: distance in Å
        :param round_digits: specify to how many decimal places round the measured distance
        :return: list of atoms at given distance from a heteroatom
        """
        atoms = []
        for atom in self.structure.get_atoms():
            if round(atom - hetatm, round_digits) == distance:
                atoms.append(atom)
        return atoms

    def get_residues_at_distance(self, hetatm, distance, round_digits=1):
        """
        Function serves to find all the residues at given distance from a heteroatom.
        :param hetatm: object representing an atom
        :param distance: distance in Å
        :param round_digits: specify to how many decimal places round the measured distance
        :return: list of residues at given distance from a heteroatom
        """
        residues = []
        for atom in self.structure.get_atoms():
            if round((hetatm - atom), round_digits) == distance:
                if atom.get_parent() not in residues:
                    residues.append(atom.get_parent())
        return residues

    def ratio_exposed(self):
        """
        Function calculates ratio of exposed and buried amino acids within a molecule using half-sphere exposure values.
        :return: string describing ratio of exposed and buried residues in %
        """
        exposed = 0
        buried = 0
        for i in(HSExposure.ExposureCN(self.structure[0])):
            if i[1] >= 27:
                 buried+=1
            else:
                exposed+=1
        return str(round(buried * 100 / (buried + exposed))) + "% of residues is buried\n"\
        + str(round(exposed * 100 / (buried + exposed))) + "% of residues is exposed"

    def get_exposed_buried(self):
        """
        Help function providing data for histogram in form of dictionaries containing numbers of each amino acids found
        in a molecule.
        :return: dictionaries for buried and for exposed amino acid compositions data
        """
        his_suf = {}
        his_bur = {}
        for i in (HSExposure.ExposureCN(self.structure[0])):
            if i[1] >= 27:
                if i[0].resname not in his_bur:
                    his_bur[i[0].resname] = 0
                his_bur[i[0].resname] += 1
            else:
                if i[0].resname not in his_suf:
                    his_suf[i[0].resname] = 0
                his_suf[i[0].resname] += 1

        return his_bur, his_suf

    def polarity_ratio(self):
        """
        Function calculates proportion of polar and nonpolar buried and exposed amino acids within a molecule.
        :return: string describing ratio of polar and nonpolar residues in %
        """
        data = self.get_exposed_buried()
        polar = 0
        nonpolar = 0
        res = ""
        for aa in data[0]:
            if aa in POLAR:
                polar += data[0][aa]
            else:
                nonpolar += data[0][aa]
        res += str(round(polar * 100 / (polar + nonpolar))) + "% of buried residues is polar\n"
        res += str(round(nonpolar * 100 / (polar + nonpolar))) + "% of buried residues is not polar\n"
        polar = 0
        nonpolar = 0
        for aa in data[1]:
            if aa in POLAR:
                polar += data[1][aa]
            else:
                nonpolar += data[1][aa]
        res += str(round(polar * 100 / (polar + nonpolar))) + "% of exposed residues is polar\n"
        res += str(round(nonpolar * 100 / (polar + nonpolar))) + "% of exposed residues is not polar"
        return res

    def show_histogram(self, exposed=False):
        """
        Function outputs a histogram showing amino acid composition of either buried or exposed residues.
        :param exposed: specifies which data to show (buried or exposed amino acid composition)
        """
        data = self.get_exposed_buried()
        i = 0
        if exposed:
            i = 1
        plt.bar(data[i].keys(), data[i].values())
        plt.show()

    def save_histogram(self, file="histogram", exposed=False):
        """
        Function stores a histogram showing amino acid composition of either buried or exposed residues.
        :param file: name of the histogram file
        :param exposed: specifies which data to show (buried or exposed amino acid composition)
        """
        data = self.get_exposed_buried()
        i = 0
        if exposed:
            i = 1
        plt.bar(data[i].keys(), data[i].values())
        plt.savefig(file)
