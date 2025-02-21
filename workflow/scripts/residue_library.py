"""
This module contains the main residue library for this repo.
"""

from rdkit import Chem

class Residue:
    def __init__(self, res_let: str, flavor: int):
        self.mol = Chem.rdmolfiles.MolFromSequence(res_let, flavor=flavor)
        self.res_name = self.mol.GetAtomWithIdx(0).GetPDBResidueInfo().GetResidueName().strip()
        # Give the OP2 atom of DNA/RNA a negative charge, removing its implicit hydrogen.
        if self.res_name in ['A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT']:
            for atom in self.mol.GetAtoms():
                if atom.GetPDBResidueInfo().GetName().strip() == 'OP2':
                    atom.SetFormalCharge(-1)
                    atom.SetNoImplicit(True)
        self.don_idx_list = _donors(self.mol)
        self.acc_idx_list = _acceptors(self.mol)


# Consider a RDKit molecule and return the indices of atoms that are hydrogen bond donors.
def _donors(mol):
    idx_list = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in [7, 8, 16] and atom.GetTotalNumHs() > 0:
            idx_list.append(atom.GetIdx())
        else:
            continue
    return idx_list


# Consider a RDKit molecule and return the indices of atoms that are hydrogen bond acceptors.
def _acceptors(mol):
    idx_list = []
    for atom in mol.GetAtoms():
        # SP3 N with 3 covalent bonds
        if atom.GetAtomicNum() == 7 and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3 and atom.GetTotalValence() == 3:
            idx_list.append(atom.GetIdx())
        # SP2 N with 2 covalent bonds
        if atom.GetAtomicNum() == 7 and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2 and atom.GetTotalValence() == 2:
            idx_list.append(atom.GetIdx())
        # oxygen or sulfur
        if atom.GetAtomicNum() in [8, 16]:
            idx_list.append(atom.GetIdx())
        else:
            continue
    return idx_list

# Define a list of the canonical protein and RNA/DNA residues, with details provided by the RDKit library.
RESIDUE_LIBRARY = [Residue('R', 0), Residue('H', 0), Residue('K', 0), Residue('D', 0), Residue('E', 0), Residue('S', 0),
                   Residue('T', 0), Residue('N', 0), Residue('Q', 0), Residue('C', 0), Residue('G', 0), Residue('P', 0),
                   Residue('A', 0), Residue('V', 0), Residue('I', 0), Residue('L', 0), Residue('M', 0), Residue('F', 0),
                   Residue('Y', 0), Residue('W', 0), Residue('A', 3), Residue('C', 3), Residue('G', 3), Residue('U', 3),
                   Residue('A', 7), Residue('C', 7), Residue('G', 7), Residue('T', 7)]

# Define a list of dictionaries that provides information on the canonical protein and RNA/DNA residues. The indices of
# each donor or acceptor atom list specifies the following:
#     0: donor/acceptor atom name
#     1: donor/acceptor element
#     2: is the donor/acceptor rotatable?
#     3: donor/acceptor formal charge according to PyMOL
#     4: antecedent atom of donor/acceptor
OLD_RESIDUE_LIBRARY = [
    {
        "res": "ALA",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"], ["OXT", "O", False, 0, "C"]]
    },
    {
        "res": "ASP",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"], ["OXT", "O", False, 0, "C"], ["OD1", "O", False, 0, "CG"],
                ["OD2", "O", False, -1, "CG"]]
    },
    {
        "res": "ASN",
        "don": [["N", "N", False, 0, "CA"], ["ND2", "N", False, 0, "CG"]],
        "acc": [["O", "O", False, 0, "C"], ["OXT", "O", False, 0, "C"], ["OD1", "O", False, 0, "CG"]]
    },
    {
        "res": "ARG",
        "don": [["N", "N", False, 0, "CA"], ["NE", "N", False, 0, "CD"], ["NH1", "N", False, 1, "CZ"],
                ["NH2", "N", False, 0, "CZ"]],
        "acc": [["O", "O", False, 0, "C"], ["OXT", "O", False, 0, "C"]]
    },
    {
        "res": "CYS",
        "don": [["N", "N", False, 0, "CA"], ["SG", "S", True, 0, "CB"]],
        "acc": [["O", "O", False, 0, "C"], ["OXT", "O", False, 0, "C"], ["SG", "S", True, 0, "CB"]]
    },
    {
        "res": "GLU",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"], ["OXT", "O", False, 0, "C"], ["OE1", "O", False, 0, "CD"],
                ["OE2", "O", False, -1, "CD"]]
    },
    {
        "res": "GLN",
        "don": [["N", "N", False, 0, "CA"], ["NE2", "N", False, 0, "CD"]],
        "acc": [["O", "O", False, 0, "C"], ["OXT", "O", False, 0, "C"], ["OE1", "O", False, 0, "CD"]]
    },
    {
        "res": "GLY",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"], ["OXT", "O", False, 0, "C"]]
    },
    {
        "res": "HIS",
        "don": [["N", "N", False, 0, "CA"], ["ND1", "N", False, 0, "CG"], ["NE2", "N", False, 0, "CD2"]],
        "acc": [["O", "O", False, 0, "C"], ["OXT", "O", False, 0, "C"], ["ND1", "N", False, 0, "CG"],
                ["NE2", "N", False, 0, "CD2"]]
    },
    {
        "res": "ILE",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"], ["OXT", "O", False, 0, "C"]]
    },
    {
        "res": "LEU",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"], ["OXT", "O", False, 0, "C"]]
    },
    {
        "res": "LYS",
        "don": [["N", "N", False, 0, "CA"], ["NZ", "N", True, 1, "CE"]],
        "acc": [["O", "O", False, 0, "C"], ["OXT", "O", False, 0, "C"]]
    },
    {
        "res": "MET",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"], ["OXT", "O", False, 0, "C"], ["SD", "S", True, 0, "CG"]]
    },
    {
        "res": "PHE",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"], ["OXT", "O", False, 0, "C"]]
    },
    {
        "res": "PRO",
        "don": [],
        "acc": [["O", "O", False, 0, "C"], ["OXT", "O", False, 0, "C"]]
    },
    {
        "res": "SER",
        "don": [["N", "N", False, 0, "CA"], ["OG", "O", True, 0, "CB"]],
        "acc": [["O", "O", False, 0, "C"], ["OXT", "O", False, 0, "C"], ["OG", "O", True, 0, "CB"]]
    },
    {
        "res": "THR",
        "don": [["N", "N", False, 0, "CA"], ["OG1", "O", True, 0, "CB"]],
        "acc": [["O", "O", False, 0, "C"], ["OXT", "O", False, 0, "C"], ["OG1", "O", True, 0, "CB"]]
    },
    {
        "res": "TRP",
        "don": [["N", "N", False, 0, "CA"], ["NE1", "N", False, 0, "CD1"]],
        "acc": [["O", "O", False, 0, "C"], ["OXT", "O", False, 0, "C"]]
    },
    {
        "res": "TYR",
        "don": [["N", "N", False, 0, "CA"], ["OH", "O", False, 0, "CZ"]],
        "acc": [["O", "O", False, 0, "C"], ["OXT", "O", False, 0, "C"], ["OH", "O", False, 0, "CZ"]]
    },
    {
        "res": "VAL",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"], ["OXT", "O", False, 0, "C"]]
    },
    {
        "res": "A",
        "don": [["O2'", "O", True, 0, "C2'"], ["N6", "N", False, 0, "C6"]],
        "acc": [["O2'", "O", True, 0, "C2'"], ["O3'", "O", False, 0, "C3'"], ["O4'", "O", False, 0, "C4'"],
                ["O5'", "O", False, 0, "C5'"], ["OP1", "O", False, 0, "P"], ["OP2", "O", False, -1, "P"],
                ["N1", "N", False, 0, "C2"], ["N3", "N", False, 0, "C2"], ["N7", "N", False, 0, "C5"]]
    },
    {
        "res": "C",
        "don": [["O2'", "O", True, 0, "C2'"], ["N4", "N", False, 0, "C4"]],
        "acc": [["O2'", "O", True, 0, "C2'"], ["O3'", "O", False, 0, "C3'"], ["O4'", "O", False, 0, "C4'"],
                ["O5'", "O", False, 0, "C5'"], ["OP1", "O", False, 0, "P"], ["OP2", "O", False, -1, "P"],
                ["O2", "O", False, 0, "C2"], ["N3", "N", False, 0, "C2"]]
    },
    {
        "res": "G",
        "don": [["O2'", "O", True, 0, "C2'"], ["N1", "N", False, 0, "C2"], ["N2", "N", False, 0, "C2"]],
        "acc": [["O2'", "O", True, 0, "C2'"], ["O3'", "O", False, 0, "C3'"], ["O4'", "O", False, 0, "C4'"],
                ["O5'", "O", False, 0, "C5'"], ["OP1", "O", False, 0, "P"], ["OP2", "O", False, -1, "P"],
                ["N3", "N", False, 0, "C2"], ["O6", "O", False, 0, "C6"], ["N7", "N", False, 0, "C5"]]
    },
    {
        "res": "U",
        "don": [["O2'", "O", True, 0, "C2'"], ["N3", "N", False, 0, "C2"]],
        "acc": [["O2'", "O", True, 0, "C2'"], ["O3'", "O", False, 0, "C3'"], ["O4'", "O", False, 0, "C4'"],
                ["O5'", "O", False, 0, "C5'"], ["OP1", "O", False, 0, "P"], ["OP2", "O", False, -1, "P"],
                ["O2", "O", False, 0, "C2"], ["O4", "O", False, 0, "C4"]]
    },
    {
        "res": "DA",
        "don": [["N6", "N", False, 0, "C6"]],
        "acc": [["O3'", "O", False, 0, "C3'"], ["O4'", "O", False, 0, "C4'"],
                ["O5'", "O", False, 0, "C5'"], ["OP1", "O", False, 0, "P"], ["OP2", "O", False, -1, "P"],
                ["N1", "N", False, 0, "C2"], ["N3", "N", False, 0, "C2"], ["N7", "N", False, 0, "C5"]]
    },
    {
        "res": "DC",
        "don": [["N4", "N", False, 0, "C4"]],
        "acc": [["O3'", "O", False, 0, "C3'"], ["O4'", "O", False, 0, "C4'"],
                ["O5'", "O", False, 0, "C5'"], ["OP1", "O", False, 0, "P"], ["OP2", "O", False, -1, "P"],
                ["O2", "O", False, 0, "C2"], ["N3", "N", False, 0, "C2"]]
    },
    {
        "res": "DG",
        "don": [["N1", "N", False, 0, "C2"], ["N2", "N", False, 0, "C2"]],
        "acc": [["O3'", "O", False, 0, "C3'"], ["O4'", "O", False, 0, "C4'"],
                ["O5'", "O", False, 0, "C5'"], ["OP1", "O", False, 0, "P"], ["OP2", "O", False, -1, "P"],
                ["N3", "N", False, 0, "C2"], ["O6", "O", False, 0, "C6"], ["N7", "N", False, 0, "C5"]]
    },
    {
        "res": "DT",
        "don": [["N3", "N", False, 0, "C2"]],
        "acc": [["O3'", "O", False, 0, "C3'"], ["O4'", "O", False, 0, "C4'"],
                ["O5'", "O", False, 0, "C5'"], ["OP1", "O", False, 0, "P"], ["OP2", "O", False, -1, "P"],
                ["O2", "O", False, 0, "C2"], ["O4", "O", False, 0, "C4"]]
    }
]
