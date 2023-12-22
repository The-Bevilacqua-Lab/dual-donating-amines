"""
This module will eventaully do something.
"""

import sys
from pymol import cmd
from pymol import stored


def eval_H_bonding(donor_index, acceptor_index, pdb):
    # initially assume that the evaluation will complete successfully
    successful_completion = True
    # initialize an empty list that can be used to document why an evaluation was not completed successfully
    notes = []
    # ensure that the indices for the donor and acceptor atoms account for exactly two atoms
    if cmd.count_atoms(f'index {donor_index} index {acceptor_index}') != 2:
        successful_completion = False
        print(f"Error: The indices provided for a potential H-bonding donor and acceptor atom pair does not account "
              f"for exactly two atoms in PDB ID {pdb}.")
        notes.append(f"Error: The indices provided for a potential H-bonding donor and acceptor atom pair does not "
                     f"account for exactly two atoms in PDB ID {pdb}.")
        return [successful_completion, notes]
    # collect information on the donor and acceptor atoms
    stored.donor_atom = ()
    cmd.iterate(f'index {donor_index}', 'stored.donor_atom = (index, name, resn, resi, chain)')
    stored.acceptor_atom = ()
    cmd.iterate(f'index {acceptor_index}', 'stored.acceptor_atom = (index, name, resn, resi, chain)')
    # construct lists containing the names of the canonical protein and nucleic residues
    protein_residues = ['ALA', 'ASP', 'ASN', 'ARG', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
                        'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    nucleic_residues = ['A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT']
    # determine whether the donor atom is at the end of a chain and is capable of donating an H-bond
    terminal_donating_atom = []
    if cmd.count_atoms(f'neighbor index {donor_index}') == 1 and (
            stored.donor_atom[1] == "N" or
            stored.donor_atom[1] == "O5'" or
            stored.donor_atom[1] == "O3'"):
        # an amino acid N atom at the end of the chain is capable of donating an H-bond
        if stored.donor_atom[1] == "N" and stored.donor_atom[2] in protein_residues:
            terminal_donating_atom = ["N", "N", True, 0, "CA"]
        else:
            stored.bonded_atom = ''
            cmd.iterate(f'neighbor index {donor_index}','stored.bonded_atom = name')
            # only an RNA O5' atom at the 5' end of the chain is capable of donating an H-bond
            if (stored.donor_atom[1] == "O5'" and stored.bonded_atom == "C5'" and stored.donor_atom[2] in
                    nucleic_residues):
                terminal_donating_atom = ["O5'", "O", True, 0, "C5'"]
            # only an RNA O3' atom at the 3' end of the chain is capable of donating an H-bond
            if (stored.donor_atom[1] == "O3'" and stored.bonded_atom == "C3'" and stored.donor_atom[2] in
                    nucleic_residues):
                terminal_donating_atom = ["O3'", "O", True, 0, "C3'"]
    # determine whether the donor atom is listed in the residue library
    # if found, collect other information about that atom
    donor_info = []
    donor_res_found = False
    donor_atom_found = False
    for residue in residue_library:
        if stored.donor_atom[2] == residue['res']:
            donor_res_found = True
            for atom in residue['don']:
                if stored.donor_atom[1] == atom[0]:
                    donor_atom_found = True
                    donor_info = atom
            if (terminal_donating_atom and (stored.donor_atom[1] == "N" and stored.donor_atom[2] in protein_residues) or
                    ((stored.donor_atom[1] == "O5'" or stored.donor_atom[1] == "O3'") and stored.donor_atom[2] in
                     nucleic_residues)):
                donor_atom_found = True
                donor_info = terminal_donating_atom
    if not donor_res_found:
        successful_completion = False
        notes.append(f"Residue {stored.donor_atom[2]} of the donor atom was not found in the residue library when "
                     f"evaluating a potential H-bond in PDB ID {pdb}.")
    if not donor_atom_found:
        successful_completion = False
        notes.append(f"The donor atom {stored.donor_atom[1]} of residue {stored.donor_atom[2]} was not found in the "
                     f"residue library when evaluating a potential H-bond in PDB ID {pdb}.")
    # determine whether the acceptor atom is listed in the residue library
    # if found, collect other information about that atom
    acceptor_info = []
    acceptor_res_found = False
    acceptor_atom_found = False
    for residue in residue_library:
        if stored.acceptor_atom[2] == residue['res']:
            acceptor_res_found = True
            for atom in residue['acc']:
                if stored.acceptor_atom[1] == atom[0]:
                    acceptor_atom_found = True
                    acceptor_info = atom
    if not acceptor_res_found:
        successful_completion = False
        notes.append(f"Residue {stored.acceptor_atom[2]} of the acceptor atom was not found in the residue library "
                     f"when evaluating a potential H-bond in PDB ID {pdb}.")
    if not acceptor_atom_found:
        successful_completion = False
        notes.append(f"The acceptor atom {stored.acceptor_atom[1]} of residue {stored.acceptor_atom[2]} was not found "
                     f"in the residue library when evaluating a potential H-bond in PDB ID {pdb}.")
    # if any of the residues of the donor or acceptor atoms or the atoms themselves were not found in the residue
    # library, return a non-successful completion with the relevant notes
    if not successful_completion:
        return [False, notes]
    # add a hydrogen to non-rotatable donors
    # if not donor_info[2] and not terminal_donating_atom:
    if not donor_info[2]:
        cmd.h_add(f'index {donor_index}')

    print(donor_info)
# region library
# define a list of dictionaries that provides information on the canonical protein and RNA/DNA residues
# the indices of each donor or acceptor atom list specifies the following:
#     0: donor/acceptor atom name
#     1: donor/acceptor element
#     2: is the donor/acceptor rotatable?
#     3: donor/acceptor formal charge
#     4: antecedent atom of donor/acceptor
residue_library = [
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
# endregion

# print(eval_H_bonding(51045, 53543, '6XU8'))
# print(eval_H_bonding(51037, 53543, '6XU8'))
# TODO look into the error message after running the following call to eval_H_bonding
print(eval_H_bonding(51059, 53543, '6XU8'))
