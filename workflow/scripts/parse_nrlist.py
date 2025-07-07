"""
This script provides a function that iterates through the lines of a representative set file and collects the
equivalence class names, PDB IDs, model info, and chain info of the equivalence class members that serve as the
representative structures of their corresponding equivalence classes. Information for the equivalence class members are
returned as a single list of strings. The name of the representative set file must be provided as the argument. Lastly,
the function also prints the number of IFEs and unique PDB IDs present in the representative set file.
"""

import sys
import csv
import numpy as np


def parse_nrlist(nrlist_file):
    # Count the number of IFEs and unique PDB IDs in the provided nrlist file.
    pdb_id_list_count = []
    with open(nrlist_file, mode='r') as read_file:
        for line in csv.reader(read_file):
            pdb_id_chars = []
            for char in line[1].split('+')[0]:
                if char != '|':
                    pdb_id_chars.append(char)
                else:
                    break
            # check if the PDB ID contains the expected number of characters
            if len("".join(pdb_id_chars)) != 4:
                print(f"Error: At least one of the PDB IDs for a member of equivalence class {line[0]} does "
                      f"not contain the expected number of characters.")
                sys.exit(1)
            pdb_id_list_count.append("".join(pdb_id_chars))
    print(f"The provided representative structure file contains {len(np.array(pdb_id_list_count))} IFEs and "
          f"{len(np.unique(np.array(pdb_id_list_count)))} unique PDB IDs.")
    # Prepare a dictionary of the equivalence classes that includes the PDB ID, model info, and chain info for the
    # equivalence class members.
    eq_class_dict = {}
    with open(nrlist_file, mode='r') as read_file:
        for line in csv.reader(read_file):
            eq_class_dict[line[0]] = {'PDB_ID': [], 'model': [], 'chain_list': []}
            pdb_id_list = []
            model_list = []
            chain_list = []
            for group in line[1].split('+'):
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
            eq_class_dict[line[0]]['PDB_ID'].append(pdb_id_list[0])
            eq_class_dict[line[0]]['model'].append(model_list[0])
            eq_class_dict[line[0]]['chain_list'].append(chain_list)
    # Prepare a list of strings where each string includes the equivalence class name, PDB ID, model info, and chain
    # info for the equivalence class members.
    eq_class_members = []
    for eq_class in eq_class_dict:
        if not (len(eq_class_dict[eq_class]['PDB_ID']) == len(eq_class_dict[eq_class]['model']) ==
                len(eq_class_dict[eq_class]['chain_list'])):
            print(f"Error: There is an issue with the information for a member of equivalence class {eq_class} in the "
                  f"Snakefile.")
            sys.exit(1)
        for idx in range(len(eq_class_dict[eq_class]['PDB_ID'])):
            eq_class_members.append(f'{eq_class}_{eq_class_dict[eq_class]["PDB_ID"][idx]}_'
                                    f'{eq_class_dict[eq_class]["model"][idx]}')
            for chain in eq_class_dict[eq_class]["chain_list"][idx]:
                eq_class_members[-1] = eq_class_members[-1] + "_" + chain
    return eq_class_members
