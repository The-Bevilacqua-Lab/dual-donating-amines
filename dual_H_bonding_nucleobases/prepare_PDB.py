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

# directory = os.getcwd()
# original_mmCIF_directory = directory + "/original_mmCIF_files"
# modified_mmCIF_directory = directory + "/modified_mmCIF_files"
#
# # identify the representative set file
# nrlist_file = parse_nrlist.identify()
#
# # iterate through the lines of the representative set file and collect the equivalence class names and PDB IDs, model
# # info, and chain info for the representative structures
# nrlist_info = parse_nrlist.get_info(nrlist_file)
#
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
donor_atoms = []
acceptor_atoms = []
for residue in eval_H_bonding.residue_library:
    for donor in residue['don']:
        donor_atoms.append((residue['res'], donor[0]))
    for acceptor in residue['acc']:
        acceptor_atoms.append((residue['res'], acceptor[0]))

# construct two strings that can be used with PyMOL to select all possible donor or all possible acceptor atoms
acceptor_string = ''
last_atom = acceptor_atoms[-1]
for atom in acceptor_atoms:
    acceptor_string += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        acceptor_string += ' '
donor_string = ''
last_atom = donor_atoms[-1]
for atom in donor_atoms:
    donor_string += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        donor_string += ' '

# # work with each equivalence class specified in the representative set file
# for eq_class in nrlist_info:
#     rep_struct_list = []
#     for i in range(len(eq_class[2])):
#         # the create function below assumes that the different representative RNA chains belong to the same model
#         # if this is not the case, print an error message and exit
#         if eq_class[2][i] != eq_class[2][0]:
#             print(f"Error: There are representative RNA chains that belong to different models for equivalence class "
#                   f"{eq_class[0]}.")
#             sys.exit(1)
#         rep_struct_list.append(f"chain {eq_class[3][i]}")
#     rep_struct = " ".join(rep_struct_list)
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
#     cmd.save(f'{modified_mmCIF_directory}/{eq_class[0]}.cif')
#     cmd.delete('all')
