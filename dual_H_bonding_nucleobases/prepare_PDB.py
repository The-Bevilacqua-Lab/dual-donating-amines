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
    stored.alt_conf = []
    cmd.iterate('not alt "" and q<1.0', 'print(ID)')





    # cmd.delete('all')
