"""
This module contains various constants that are used by other modules in this package.
"""

# construct three tuples of tuples containing atom and residue names that describe either donor, protonated donor
# (formal charge of +1), or acceptor atoms of particular interest
DONORS_OF_INTEREST = (('A', 'N6'), ('C', 'N4'), ('G', 'N2'))
PROT_DONORS_OF_INTEREST = (('A', 'N1'), ('A', 'N3'), ('A', 'N7'), ('C', 'N3'), ('G', 'N3'), ('G', 'N7'))
ACCEPTORS_OF_INTEREST = (('A', 'N1'), ('A', 'N3'), ('A', 'N7'), ('C', 'O2'), ('C', 'N3'), ('G', 'N3'), ('G', 'O6'),
                         ('G', 'N7'))
DEPROT_ACCEPTORS_OF_INTEREST = (('G', 'N1'), ('U', 'N3'))

# define a list of dictionaries that provides information on the canonical protein and RNA/DNA residues
# the indices of each donor or acceptor atom list specifies the following:
#     0: donor/acceptor atom name
#     1: donor/acceptor element
#     2: is the donor/acceptor rotatable?
#     3: donor/acceptor formal charge
#     4: antecedent atom of donor/acceptor
RESIDUE_LIBRARY = [
    {
        "res": "ALA",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"]]
    },
    {
        "res": "ASP",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"], ["OD1", "O", False, 0, "CG"], ["OD2", "O", False, -1, "CG"]]
    },
    {
        "res": "ASN",
        "don": [["N", "N", False, 0, "CA"], ["ND2", "N", False, 0, "CG"]],
        "acc": [["O", "O", False, 0, "C"], ["OD1", "O", False, 0, "CG"]]
    },
    {
        "res": "ARG",
        "don": [["N", "N", False, 0, "CA"], ["NE", "N", False, 0, "CD"], ["NH1", "N", False, 1, "CZ"],
                ["NH2", "N", False, 0, "CZ"]],
        "acc": [["O", "O", False, 0, "C"]]
    },
    {
        "res": "CYS",
        "don": [["N", "N", False, 0, "CA"], ["SG", "S", True, 0, "CB"]],
        "acc": [["O", "O", False, 0, "C"], ["SG", "S", True, 0, "CB"]]
    },
    {
        "res": "GLU",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"], ["OE1", "O", False, 0, "CD"], ["OE2", "O", False, -1, "CD"]]
    },
    {
        "res": "GLN",
        "don": [["N", "N", False, 0, "CA"], ["NE2", "N", False, 0, "CD"]],
        "acc": [["O", "O", False, 0, "C"], ["OE1", "O", False, 0, "CD"]]
    },
    {
        "res": "GLY",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"]]
    },
    {
        "res": "HIS",
        "don": [["N", "N", False, 0, "CA"], ["ND1", "N", False, 0, "CG"], ["NE2", "N", False, 0, "CD2"]],
        "acc": [["O", "O", False, 0, "C"], ["ND1", "N", False, 0, "CG"], ["NE2", "N", False, 0, "CD2"]]
    },
    {
        "res": "ILE",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"]]
    },
    {
        "res": "LEU",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"]]
    },
    {
        "res": "LYS",
        "don": [["N", "N", False, 0, "CA"], ["NZ", "N", True, 1, "CE"]],
        "acc": [["O", "O", False, 0, "C"]]
    },
    {
        "res": "MET",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"], ["SD", "S", True, 0, "CG"]]
    },
    {
        "res": "PHE",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"]]
    },
    {
        "res": "PRO",
        "don": [],
        "acc": [["O", "O", False, 0, "C"]]
    },
    {
        "res": "SER",
        "don": [["N", "N", False, 0, "CA"], ["OG", "O", True, 0, "CB"]],
        "acc": [["O", "O", False, 0, "C"], ["OG", "O", True, 0, "CB"]]
    },
    {
        "res": "THR",
        "don": [["N", "N", False, 0, "CA"], ["OG1", "O", True, 0, "CB"]],
        "acc": [["O", "O", False, 0, "C"], ["OG1", "O", True, 0, "CB"]]
    },
    {
        "res": "TRP",
        "don": [["N", "N", False, 0, "CA"], ["NE1", "N", False, 0, "CD1"]],
        "acc": [["O", "O", False, 0, "C"]]
    },
    {
        "res": "TYR",
        "don": [["N", "N", False, 0, "CA"], ["OH", "O", False, 0, "CZ"]],
        "acc": [["O", "O", False, 0, "C"], ["OH", "O", False, 0, "CZ"]]
    },
    {
        "res": "VAL",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"]]
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
