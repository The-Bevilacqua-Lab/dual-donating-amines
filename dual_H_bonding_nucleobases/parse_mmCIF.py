"""
This module reads the representative set file and iterates through the mmCIF files to collect information related to
each equivalence class that may be useful for data analysis.
"""

import sys
import os
import csv

directory = os.getcwd()

# initialize a dictionary that will store information from the representative set file and mmCIF files
eq_class_info = {'equivalence_class': [], 'nrlist_PDB_ID': [], 'model': [], 'chain': [], 'file_name': [],
                 'mmCIF_PDB_ID': [], 'release_date': [], 'method': [], 'resolution': [], 'organism': []}

# identify the representative set file
nrlist_file = ""
try:
    file_list = os.listdir(directory)
    num_nrlist = 0
    for file in file_list:
        if file[0:6] == "nrlist":
            if file[-3:] == "csv":
                num_nrlist += 1
                nrlist_file = file
    if num_nrlist == 0:
        raise Exception("There is no representative set file in the current working directory.")
    if num_nrlist > 1:
        raise Exception("There is more than one representative set file in the current working directory.")
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)

with open(nrlist_file, mode='r') as file:
    for line in csv.reader(file):
        PDB_ID_list = []
        model_list = []
        chain_list = []
        for group in line[1].split('+'):
            PDB_ID = []
            model = []
            chain = []
            index = 0
            for char in group:
                if char != '|' and index == 0:
                    PDB_ID.append(char)
                elif char != '|' and index == 1:
                    model.append(char)
                elif char != '|' and index == 2:
                    chain.append(char)
                elif char == '|':
                    index += 1
            # check if the PDB ID contains the expected number of characters
            if len("".join(PDB_ID)) != 4:
                print("Error: At least one of the PDB IDs does not contain the expected number of characters.")
                sys.exit(1)
            PDB_ID_list.append("".join(PDB_ID))
            model_list.append("".join(model))
            chain_list.append("".join(chain))
        eq_class_info['equivalence_class'].append(line[0])
        eq_class_info['nrlist_PDB_ID'].append(PDB_ID_list)
        eq_class_info['model'].append(model_list)
        eq_class_info['chain'].append(chain_list)


# mmCIF_files = []
#
# try:
#     mmCIF_files = os.listdir(directory + "/mmCIF_files")
# except FileNotFoundError:
#     print("Error: There is no folder named mmCIF_files in the current working directory.")
#     sys.exit(1)
# except Exception as e:
#     print(f"Error: {e}")
#     sys.exit(1)



# for f in mmCIF_files:
#     if f[-4:] != ".cif":
#         continue
#     else:

