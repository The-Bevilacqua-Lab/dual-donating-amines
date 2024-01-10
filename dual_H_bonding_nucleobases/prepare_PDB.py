"""
This module reads information on the equivalence classes in the representative set file. For each equivalence class, it
works with PyMOL to save a mmCIF file of the entire structure containing the representative RNA chain(s) to the
mmCIF_files folder and a PDB file of either the entire structure or a portion of the structure containing the
representative RNA chain(s) to the PDB_files folder. In addition to the representative RNA chain(s), the PDB file
contains other chains that are close to the representative RNA chain(s). Also, alternate conformations are removed such
that each atom in the resulting PDB file only has a single conformation.
"""

import sys
import os
from pymol import cmd
from pymol import stored
import parse_nrlist

directory = os.getcwd()
original_mmCIF_directory = directory + "/original_mmCIF_files"
modified_mmCIF_directory = directory + "/modified_mmCIF_files"

# TODO atom IDs may not be unique, modify the code to insure the atom selections are specific enough
# identify the representative set file
nrlist_file = parse_nrlist.identify()

# iterate through the lines of the representative set file and collect the equivalence class names and PDB IDs, model
# info, and chain info for the representative structures
nrlist_info = parse_nrlist.get_info(nrlist_file)

# check whether an original_mmCIF_files folder exists
# if it exists and contains files, exit with an error message
# if it does not exist, create the folder
if os.path.isdir(original_mmCIF_directory):
    if len(os.listdir(original_mmCIF_directory)) != 0:
        print("Error: There are already files in the original_mmCIF_files folder. This folder should be empty prior to "
              "running the prepare_PDB module.")
        sys.exit(1)
else:
    os.mkdir(original_mmCIF_directory)

# change fetch_path to the original_mmCIF_files folder in the current working directory so that PyMOL drops the mmCIF
# files into this folder when running fetch
cmd.set('fetch_path', cmd.exp_path('./original_mmCIF_files'))

# check whether a modified_mmCIF_files folder exists
# if it exists and contains files, exit with an error message
# if it does not exist, create the folder
if os.path.isdir(modified_mmCIF_directory):
    if len(os.listdir(modified_mmCIF_directory)) != 0:
        print("Error: There are already files in the modified_mmCIF_files folder. This folder should be empty prior to "
              "running the prepare_PDB module.")
        sys.exit(1)
else:
    os.mkdir(modified_mmCIF_directory)

# work with each equivalence class specified in the representative set file
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
    # delete all objects in the current PyMOL session
    cmd.delete('all')
    # retrieve the structure that contains the representative RNA chains
    cmd.fetch(eq_class[1][0])
    # if the structure contains more than one model (or "state" in PyMOL), create a new PyMOL object of just that model
    if cmd.count_states(eq_class[1][0]) > 1:
        cmd.create(f'{eq_class[1][0]}_state_{eq_class[2][0]}', selection=eq_class[1][0], source_state=eq_class[2][0],
                   target_state=1)
        cmd.delete(eq_class[1][0])
    # store atoms that have alternate conformations
    stored.alt_atoms = []
    cmd.iterate_state(0, 'not alt ""', 'stored.alt_atoms.append((name,resn,resi,chain,alt,q))')
    # prepare a list of atom groups
    # each atom group consists of one or more different conformations of the same atom
    alt_atoms_grouped = []
    for atom in stored.alt_atoms:
        match = False
        for group in alt_atoms_grouped:
            if atom[0] == group[0][0] and atom[1] == group[0][1] and atom[2] == group[0][2] and atom[3] == group[0][3]:
                group.append(atom)
                match = True
        if not match:
            alt_atoms_grouped.append([atom])
    # determine which atom conformations to keep, prioritizing highest occupancy factor and lowest alt ID, in that order
    atoms_to_remove = []
    for group in alt_atoms_grouped:
        keep = []
        remove = []
        for atom in group:
            if len(keep) == 0:
                keep = atom
            elif atom[5] > keep[5]:
                remove.append(keep)
                keep = atom
            elif atom[5] == keep[5]:
                if len(atom[4]) != 1 or len(keep[4]) != 1:
                    print(f"Error: There is at least one atom in PDB ID {eq_class[1][0]} with an alt ID that does not "
                          f"consist of exactly one character.")
                    sys.exit(1)
                elif atom[4] < keep[4]:
                    remove.append(keep)
                    keep = atom
                elif atom[4] == keep[4]:
                    print(f"Error: There are multiple atoms with the same alt ID that should have different alt IDs in "
                          f"PDB ID {eq_class[1][0]}.")
                    sys.exit(1)
                elif atom[4] > keep[4]:
                    remove.append(atom)
            elif atom[5] < keep[5]:
                remove.append(atom)
        for atom in remove:
            atoms_to_remove.append(atom)
    # TODO add code to test that bonded atoms either do not have an alt conformation or are of the same alt conformation
    # remove the atom conformations that are not going to be kept
    for atom in atoms_to_remove:
        if cmd.count_atoms(f'name {atom[0]} and resn {atom[1]} and resi {atom[2]} and chain {atom[3]} and '
                           f'alt {atom[4]}') != 1:
            print(f"Error: The info provided for an atom to remove does not account for exactly one atom in PDB ID "
                  f"{eq_class[1][0]}.")
            sys.exit(1)
        else:
            cmd.remove(f'name {atom[0]} and resn {atom[1]} and resi {atom[2]} and chain {atom[3]} and alt {atom[4]}')
    cmd.save(f'{modified_mmCIF_directory}/{eq_class[0]}.cif')
    cmd.delete('all')
