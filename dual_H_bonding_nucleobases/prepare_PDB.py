"""
This module will eventually do something.
"""

import sys
import os
import csv
from pymol import cmd
from pymol import stored
import parse_nrlist

directory = os.getcwd()
mmCIF_directory = directory + "/mmCIF_files"
PDB_directory = directory + "/PDB_files"

# identify the representative set file
nrlist_file = parse_nrlist.identify()

# iterate through the lines of the representative set file and collect the equivalence class names and PDB IDs, model
# info, and chain info for the representative structures
nrlist_info = parse_nrlist.get_info(nrlist_file)

# check whether a mmCIF_files folder exists
# if it exists and contains files, exit with an error message
# if it does not exist, create the folder
if os.path.isdir(mmCIF_directory):
    if len(os.listdir(mmCIF_directory)) != 0:
        print("Error: There are already files in the mmCIF_files folder. This folder should be empty prior to running "
              "the prepare_PDB module.")
        sys.exit(1)
else:
    os.mkdir(mmCIF_directory)

os.chdir(mmCIF_directory)

# check whether a PDB_files folder exists
# if it exists and contains files, exit with an error message
# if it does not exist, create the folder
if os.path.isdir(PDB_directory):
    if len(os.listdir(PDB_directory)) != 0:
        print("Error: There are already files in the PDB_files folder. This folder should be empty prior to running "
              "the prepare_PDB module.")
        sys.exit(1)
else:
    os.mkdir(PDB_directory)

for eq_class in nrlist_info:
    rep_struct_list = []
    for i in range(len(eq_class[2])):
        # the create function below assumes that the different representative RNA chains belong to the same model
        # if this is not the case, print an error message and exit
        if eq_class[2][i] != eq_class[2][0]:
            print(f"Error: There are representative RNA chains that belong to different models for equivalence class "
                  f"{eq_class[0]}.")
            sys.exit(1)
        rep_struct_list.append(f"chain {eq_class[3][i]}")
    rep_struct = " ".join(rep_struct_list)
    cmd.fetch(eq_class[1][0])
    # create a copy of the structure of just the model (state in PyMOL) indicated in the representative set file
    cmd.create('copy', selection=eq_class[1][0], source_state=eq_class[2][0], target_state=1)
    cmd.delete(eq_class[1][0])
    stored.alt_atoms = []
    cmd.iterate_state(0, 'not alt ""', 'stored.alt_atoms.append([ID,q,alt,state,chain,resi,name])')
    alt_atoms_grouped = []
    for atom in stored.alt_atoms:
        match = False
        for group in alt_atoms_grouped:
            if atom[3] == group[0][3] and atom[4] == group[0][4] and atom[5] == group[0][5] and atom[6] == group[0][6]:
                group.append(atom)
                match = True
        if not match:
            alt_atoms_grouped.append([atom])
    # determine which atoms to keep, prioritizing highest occupancy factor and lowest alt ID, in that order
    atoms_to_remove = []
    for group in alt_atoms_grouped:
        keep = []
        remove = []
        for atom in group:
            if len(keep) == 0:
                keep = atom
            elif atom[1] > keep[1]:
                remove.append(keep)
                keep = atom
            elif atom[1] == keep[1]:
                if len(atom[2]) != 1 or len(keep[2]) != 1:
                    print(f"Error: There is at least one atom in PDB ID {eq_class[1][0]} with an alt ID that consists "
                          f"of more than one character.")
                    sys.exit(1)
                elif atom[2] < keep[2]:
                    remove.append(keep)
                    keep = atom
                elif atom[2] == keep[2]:
                    print(f"Error: There are multiple atoms with the same alt ID that should have different alt IDs in "
                          f"PDB ID {eq_class[1][0]}.")
                    sys.exit(1)
                elif atom[2] > keep[2]:
                    remove.append(atom)
            elif atom[1] < keep[1]:
                remove.append(atom)
        for atom in remove:
            atoms_to_remove.append(atom[0])
    for atom_ID in atoms_to_remove:
        cmd.remove(f'ID {atom_ID}')
    # cmd.save(f'{PDB_directory}/{eq_class[0]}.pdb', state=eq_class[2][0])
    # cmd.delete('all')
