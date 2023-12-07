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

for eq_class in nrlist_info:
    rep_struct_list = []
    for i in range(len(eq_class[2])):
        rep_struct_list.append(f"(state {eq_class[2][i]} and chain {eq_class[3][i]})")
    rep_struct = " ".join(rep_struct_list)
    cmd.fetch(eq_class[1][0])
    cmd.remove(f'not bychain all within 5 of {rep_struct}')
    stored.alt_atoms = []
    cmd.iterate_state(0, 'not alt ""', 'stored.alt_atoms.append([ID,q,b,state,chain,resi,name])')
    alt_atoms_grouped = []
    for atom in stored.alt_atoms:
        match = False
        for group in alt_atoms_grouped:
            if atom[3] == group[0][3] and atom[4] == group[0][4] and atom[5] == group[0][5] and atom[6] == group[0][6]:
                group.append(atom)
                match = True
        if not match:
            alt_atoms_grouped.append([atom])
    # determine which atoms to keep, prioritizing highest occupancy factor, lowest b-factor, and lowest atom ID, in
    # that order
    atom_to_remove = []
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
                if atom[2] < keep[2]:
                    remove.append(keep)
                    keep = atom
                elif atom[2] == keep[2]:
                    if atom[0] < keep[0]:
                        remove.append(keep)
                        keep = atom
                    elif atom[0] == keep[0]:
                        print(f"Error: There are multiple atoms with the same ID in PDB ID {eq_class[1][0]}.")
                        sys.exit(1)
                    elif atom[0] > keep[0]:
                        remove.append(atom)
                elif atom[2] > keep[2]:
                    remove.append(atom)
            elif atom[1] < keep[1]:
                remove.append(atom)
        for atom in remove:
            atom_to_remove.append(atom[0])




    # cmd.delete('all')
