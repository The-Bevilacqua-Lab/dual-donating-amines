"""
This script provides a function that reads a full integrated functional elements (IFEs) file and collects the
equivalence class names, PDB IDs, model info, and chain info of the representative IFEs that serve as the
representative structures of their corresponding equivalence classes. Information for the representative IFEs are
returned as a single list of strings. The name of the IFEs file must be provided as the argument. Lastly, the function
also prints the number of representative IFEs and the number of unique PDB IDs among those IFEs.
"""

import sys
import csv
import numpy as np
import pandas as pd


def parse_ifes(ifes_file):
    ifes_df = pd.read_csv(ifes_file)
    # Prepare a dictionary of the equivalence classes that includes the PDB ID, model info, and chain info for the
    # representative member of each equivalence class.
    eq_class_dict = {}
    pdb_id_grand_list = []
    for idx, row in enumerate(ifes_df[ifes_df['ec_rank'] == 1].itertuples()):
        pdb_id_list = []
        model_list = []
        chain_list = []
        for group in row.ife_id.split('+'):
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
                print(f"Error: A member of equivalence class {row.ec_id} is missing information.")
                sys.exit(1)
            # check if the PDB ID contains the expected number of characters
            if len("".join(pdb_id)) != 4:
                print(f"Error: At least one of the PDB IDs for the representative member of equivalence class "
                      f"{row.ec_id} does not contain the expected number of characters.")
                sys.exit(1)
            pdb_id_list.append("".join(pdb_id))
            model_list.append("".join(model))
            chain_list.append("".join(chain))
        # Check to make sure the PDB IDs match.
        for pdb in pdb_id_list:
            if pdb != pdb_id_list[0]:
                print(f"Error: The PDB IDs for the representative member of equivalence class {row.ec_id} do not "
                      f"match.")
                sys.exit(1)
        # Check to make sure the models match.
        for model in model_list:
            if model != model_list[0]:
                print(f"Error: The models for the representative member of equivalence class {row.ec_id} do not match.")
                sys.exit(1)
        eq_class_dict[row.ec_id] = {'PDB_ID': pdb_id_list[0], 'model': model_list[0], 'chain_list': chain_list}
        pdb_id_grand_list.append(pdb_id_list[0])
    # Report the number of representative IFEs and unique PDB IDs for those IFEs in the provided IFEs file.
    print(f"The provided representative structure file contains {len(np.array(pdb_id_grand_list))} IFEs and "
          f"{len(np.unique(np.array(pdb_id_grand_list)))} unique PDB IDs.")
    # Write a filtered version of the IFEs file with only certain columns and just the representative IFEs.
    columns = ['ec_id', 'ife_id', 'pdb_resolution', 'ec_cqs', 'nts_observed', 'rfam_max_nts',
               'rfam_fraction_observed', 'rfam', 'standardized_name', 'pdb_experimental_technique']
    ifes_df.loc[(ifes_df['ec_rank'] == 1), columns].to_csv(f"{ifes_file[:-4]}_filtered.csv",
                                                           index=False, na_rep='NA')
    # Prepare a list of strings where each string includes the equivalence class name, PDB ID, model info, and chain
    # info for the representative member of the corresponding equivalence class.
    eq_class_members = []
    for eq_class in eq_class_dict:
        eq_class_members.append(f'{eq_class}_{eq_class_dict[eq_class]["PDB_ID"]}_{eq_class_dict[eq_class]["model"]}')
        for chain in eq_class_dict[eq_class]["chain_list"]:
            eq_class_members[-1] = eq_class_members[-1] + "_" + chain
    return eq_class_members
