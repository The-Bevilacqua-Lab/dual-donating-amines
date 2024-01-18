"""
This module reads information on the equivalence classes in the representative set file. For each equivalence class, it
works with PyMOL to save a mmCIF file of the original structure to the original_mmCIF_files folder and another mmCIF
file with modifications made (alternative conformations removed and hydrogens added) to the modified_mmCIF_files folder.
"""

import sys
import os
from pymol import cmd
from pymol import stored
import parse_nrlist
import eval_H_bonding
import remove_alt_conf

directory = os.getcwd()
original_mmCIF_directory = directory + "/original_mmCIF_files"
modified_mmCIF_directory = directory + "/modified_mmCIF_files"

# identify the representative set file
nrlist_file = parse_nrlist.identify()

# iterate through the lines of the representative set file and collect the equivalence class names and PDB IDs, model
# info, and chain info for the representative structures
nrlist_info = parse_nrlist.get_info(nrlist_file)

# # check whether an original_mmCIF_files folder exists
# # if it exists and contains files, exit with an error message
# # if it does not exist, create the folder
# if os.path.isdir(original_mmCIF_directory):
#     if len(os.listdir(original_mmCIF_directory)) != 0:
#         print("Error: There are already files in the original_mmCIF_files folder. This folder should be empty prior to "
#               "running the prepare_PDB module.")
#         sys.exit(1)
# else:
#     os.mkdir(original_mmCIF_directory)
#
# # change fetch_path to the original_mmCIF_files folder in the current working directory so that PyMOL drops the mmCIF
# # files into this folder when running fetch
# cmd.set('fetch_path', cmd.exp_path('./original_mmCIF_files'))
#
# # check whether a modified_mmCIF_files folder exists
# # if it exists and contains files, exit with an error message
# # if it does not exist, create the folder
# if os.path.isdir(modified_mmCIF_directory):
#     if len(os.listdir(modified_mmCIF_directory)) != 0:
#         print("Error: There are already files in the modified_mmCIF_files folder. This folder should be empty prior to "
#               "running the prepare_PDB module.")
#         sys.exit(1)
# else:
#     os.mkdir(modified_mmCIF_directory)

# construct two lists of tuples containing atom and residue names that describe either donor or acceptor atoms
# additionally, construct a list of tuples containing atom and residue names that describe rotatable donor atoms
donor_atoms = []
rotatable_donor_atoms = []
acceptor_atoms = []
for residue in eval_H_bonding.residue_library:
    for donor in residue['don']:
        donor_atoms.append((residue['res'], donor[0]))
        if donor[2]:
            rotatable_donor_atoms.append((residue['res'], donor[0]))
    for acceptor in residue['acc']:
        acceptor_atoms.append((residue['res'], acceptor[0]))

# construct three strings that can be used with PyMOL to select all possible donor, all possible rotatable donor, or
# all possible acceptor atoms
donor_string = ''
last_atom = donor_atoms[-1]
for atom in donor_atoms:
    donor_string += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        donor_string += ' '
rotatable_donor_string = ''
last_atom = rotatable_donor_atoms[-1]
for atom in rotatable_donor_atoms:
    rotatable_donor_string += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        rotatable_donor_string += ' '
acceptor_string = ''
last_atom = acceptor_atoms[-1]
for atom in acceptor_atoms:
    acceptor_string += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        acceptor_string += ' '

# construct two lists of tuples containing atom and residue names that describe either donor or acceptor atoms that
# belong to side chains that have conformations that may have been ambiguously assigned
ambiguous_donor_atoms = []
ambiguous_acceptor_atoms = []
for residue in eval_H_bonding.residue_library:
    if residue['res'] in ["ASN", "GLN", "HIS"]:
        for donor in residue['don']:
            if donor[0] != "N":
                ambiguous_donor_atoms.append((residue['res'], donor[0]))
        for acceptor in residue['acc']:
            if acceptor[0] != "O":
                ambiguous_acceptor_atoms.append((residue['res'], acceptor[0]))

# construct two strings that can be used with PyMOL to select all possible donor or all possible acceptor atoms that
# belong to side chains that have conformations that may have been ambiguously assigned
ambiguous_donor_string = ''
last_atom = ambiguous_donor_atoms[-1]
for atom in ambiguous_donor_atoms:
    ambiguous_donor_string += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        ambiguous_donor_string += ' '
ambiguous_acceptor_string = ''
last_atom = ambiguous_acceptor_atoms[-1]
for atom in ambiguous_acceptor_atoms:
    ambiguous_acceptor_string += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        ambiguous_acceptor_string += ' '

# construct two tuples of tuples containing atom and residue names that describe either donor or acceptor atoms of
# particular interest
donors_of_interest = (('A', 'N6'), ('C', 'N4'), ('G', 'N2'))
acceptors_of_interest = (('A', 'N1'), ('A', 'N3'), ('C', 'O2'), ('C', 'N3'), ('G', 'N3'), ('G', 'O6'))

# construct two strings that can be used with PyMOL to select all possible donor or all possible acceptor atoms of
# particular interest
donors_of_interest_str = ''
last_atom = donors_of_interest[-1]
for atom in donors_of_interest:
    donors_of_interest_str += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        donors_of_interest_str += ' '
acceptors_of_interest_str = ''
last_atom = acceptors_of_interest[-1]
for atom in acceptors_of_interest:
    acceptors_of_interest_str += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        acceptors_of_interest_str += ' '

# define the donor and acceptor atoms of the RNA nucleobases
adenine_nuc_donors = ["N6"]
adenine_nuc_acceptors = ["N1", "N3", "N7"]
cytosine_nuc_donors = ["N4"]
cytosine_nuc_acceptors = ["O2", "N3"]
guanine_nuc_donors = ["N1", "N2"]
guanine_nuc_acceptors = ["N3", "O6", "N7"]
uracil_nuc_donors = ["N3"]
uracil_nuc_acceptors = ["O2", "O4"]

# define the distance used to search for nearby donor or acceptor atoms
search_dist = 3.6
# define the distance used to search for nearby donor or acceptor atoms that belong to side chains that have
# conformations that may have been ambiguously assigned
search_dist_amb = search_dist + 2.5

# work with each equivalence class specified in the representative set file
for eq_class in nrlist_info:
    rep_rna_list = []
    for i in range(len(eq_class[2])):
        # the create function below assumes that the different representative RNA chains belong to the same model
        # if this is not the case, print an error message and exit
        if eq_class[2][i] != eq_class[2][0]:
            print(f"Error: There are representative RNA chains that belong to different models for equivalence class "
                  f"{eq_class[0]}.")
            sys.exit(1)
        rep_rna_list.append(f"chain {eq_class[3][i]}")
    rep_rna = " ".join(rep_rna_list)
#     # delete all objects in the current PyMOL session
#     cmd.delete('all')
#     # retrieve the structure that contains the representative RNA chains
#     cmd.fetch(eq_class[1][0])
#     # if the structure contains more than one model (or "state" in PyMOL), create a new PyMOL object of just that model
#     if cmd.count_states(eq_class[1][0]) > 1:
#         cmd.create(f'{eq_class[1][0]}_state_{eq_class[2][0]}', selection=eq_class[1][0], source_state=eq_class[2][0],
#                    target_state=1)
#         cmd.delete(eq_class[1][0])
#     # remove atoms representing alternative conformations
#     remove_status = remove_alt_conf.remove(eq_class[1][0])
#     successful_completion = remove_status[0]
#     if not successful_completion:
#         for note in remove_status[1]:
#             print(note)
#         sys.exit(1)
    # store a list of donors of interest from the representative RNA
    stored.donor_list = []
    cmd.iterate(f'{rep_rna} and ({donors_of_interest_str})',
                'stored.donor_list.append((index, name, resn, resi, chain))')
    # store a list of acceptors of interest from the representative RNA
    stored.acceptor_list = []
    cmd.iterate(f'{rep_rna} and ({acceptors_of_interest_str})',
                'stored.acceptor_list.append((index, name, resn, resi, chain))')
    # store a list of atoms near the donors of interest
    atoms_near_donors = []
    stored.donor_list = stored.donor_list[549:580]
    for donor in stored.donor_list:
        stored.nearby_atoms = []
        cmd.iterate(f'index {donor[0]} around {search_dist}',
                    'stored.nearby_atoms.append((index, name, resn, resi, chain))')
        atoms_near_donors.append(stored.nearby_atoms)
    # extract the atoms that can act as H-bond acceptors
    # they should not belong to the nucleobase of the donor
    acceptors_near_donors = []
    for atom_group in enumerate(atoms_near_donors):
        list_of_acceptors = []
        for atom in atom_group[1]:
            for acceptor in acceptor_atoms:
                if ((atom[1] == acceptor[1] and atom[2] == acceptor[0]) and not (
                        (stored.donor_list[atom_group[0]][2] == atom[2] == "A" and
                            stored.donor_list[atom_group[0]][3] == atom[3] and
                            stored.donor_list[atom_group[0]][4] == atom[4] and
                            atom[1] in adenine_nuc_acceptors) or
                        (stored.donor_list[atom_group[0]][2] == atom[2] == "C" and
                            stored.donor_list[atom_group[0]][3] == atom[3] and
                            stored.donor_list[atom_group[0]][4] == atom[4] and
                            atom[1] in cytosine_nuc_acceptors) or
                        (stored.donor_list[atom_group[0]][2] == atom[2] == "G" and
                            stored.donor_list[atom_group[0]][3] == atom[3] and
                            stored.donor_list[atom_group[0]][4] == atom[4] and
                            atom[1] in guanine_nuc_acceptors) or
                        (stored.donor_list[atom_group[0]][2] == atom[2] == "U" and
                            stored.donor_list[atom_group[0]][3] == atom[3] and
                            stored.donor_list[atom_group[0]][4] == atom[4] and
                            atom[1] in uracil_nuc_acceptors)
                        )):
                    list_of_acceptors.append(atom)
        acceptors_near_donors.append(list_of_acceptors)

    donor_h_bonds = []
    for donor in enumerate(stored.donor_list[:1]):
        h_bonds = []
        for acceptor in acceptors_near_donors[donor[0]]:
            h_bonds.append(eval_H_bonding.evaluate(donor[1], acceptor, '6xu8'))
        donor_h_bonds.append(h_bonds)

    # for x in enumerate(acceptors_near_donors):
    #     for y in x[1]:
    #         print(str(x[0]) + "\t" + str(y))
#     cmd.save(f'{modified_mmCIF_directory}/{eq_class[0]}.cif')
#     cmd.delete('all')
