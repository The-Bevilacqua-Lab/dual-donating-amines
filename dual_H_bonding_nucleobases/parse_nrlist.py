"""
This module provides functions that work with the representative set file.
"""

import sys
import os
import csv


def identify():
    """
    This function searches for and returns the name of a representative set file in the current working directory.
    If there is not exactly one representative set file, an error will be reported and the code will exit.
    """
    directory = os.getcwd()
    nrlist_file = ""
    file_list = os.listdir(directory)
    num_nrlist = 0
    for file in file_list:
        if file[0:6] == "nrlist":
            if file[-3:] == "csv":
                num_nrlist += 1
                nrlist_file = file
    if num_nrlist == 0:
        print("Error: There is no representative set file in the current working directory.")
        sys.exit(1)
    if num_nrlist > 1:
        print("Error: There is more than one representative set file in the current working directory.")
        sys.exit(1)
    return nrlist_file


def get_info(nrlist_file):
    """
    This function iterates through the lines of a representative set file and returns a list with the equivalence class
    names and PDB IDs, model info, and chain info for the representative structures. The name of the representative set
    file is a required argument.
    """
    nrlist_info = []
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
                    print("Error: At least one of the PDB IDs from the representative set file does not contain the "
                          "expected number of characters.")
                    sys.exit(1)
                PDB_ID_list.append("".join(PDB_ID))
                model_list.append("".join(model))
                chain_list.append("".join(chain))
            # check to make sure the PDB IDs match
            for pdb in PDB_ID_list:
                if pdb != PDB_ID_list[0]:
                    print("Error: The PDB IDs for the representative structure of a equivalence class do not match.")
                    sys.exit(1)
            nrlist_info.append([line[0], PDB_ID_list, model_list, chain_list])
    return nrlist_info
