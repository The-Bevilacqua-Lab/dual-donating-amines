"""
This script iterates through the lines of a representative set file and collects the equivalence class name, PDB ID,
model info, and chain info of the representative structures for each equivalence class. The info for all equivalence
classes is returned as a single dictionary. The name of the representative set file must be provided as the first
argument. If True is provided as the second argument, information for all members (not just the representative
structures) of the equivalence classes will be compiled.
"""

import sys
import csv


def parse_nrlist(nrlist_file, all_members=False):
    eq_class = {}
    with open(nrlist_file, mode='r') as read_file:
        for line in csv.reader(read_file):
            eq_class[line[0]] = {'PDB_ID': [], 'model': [], 'chain_list': []}
            for class_member in line[2].split(','):
                pdb_id_list = []
                model_list = []
                chain_list = []
                for group in class_member.split('+'):
                    pdb_id = []
                    model = []
                    chain = []
                    index = 0
                    for char in group:
                        if char != '|' and index == 0:
                            pdb_id.append(char)
                        elif char != '|' and index == 1:
                            model.append(char)
                        elif char != '|' and index == 2:
                            chain.append(char)
                        elif char == '|':
                            index += 1
                    # check to ensure this part of the representative structure contains the expected amount of
                    # information
                    if index != 2:
                        print(f"Error: A member of equivalence class {line[0]} is missing information.")
                        sys.exit(1)
                    # check if the PDB ID contains the expected number of characters
                    if len("".join(pdb_id)) != 4:
                        print(f"Error: At least one of the PDB IDs for a member of equivalence class {line[0]} does "
                              f"not contain the expected number of characters.")
                        sys.exit(1)
                    pdb_id_list.append("".join(pdb_id))
                    model_list.append("".join(model))
                    chain_list.append("".join(chain))
                # Check to make sure the PDB IDs match.
                for pdb in pdb_id_list:
                    if pdb != pdb_id_list[0]:
                        print(f"Error: The PDB IDs for a member of equivalence class {line[0]} do not match.")
                        sys.exit(1)
                # Check to make sure the models match.
                for model in model_list:
                    if model != model_list[0]:
                        print(f"Error: The models for a member of equivalence class {line[0]} do not match.")
                        sys.exit(1)
                eq_class[line[0]]['PDB_ID'].append(pdb_id_list[0])
                eq_class[line[0]]['model'].append(model_list[0])
                eq_class[line[0]]['chain_list'].append(chain_list)
                if not all_members:
                    break
    return eq_class
