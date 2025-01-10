"""
This module reads the representative set file and iterates through the mmCIF files to collect information related to
each equivalence class that may be useful for data analysis.
"""

import sys
import os
import csv
from gemmi import cif
import parse_nrlist

directory = os.getcwd()

# identify the representative set file
nrlist_file = parse_nrlist.identify()

# iterate through the lines of the representative set file and collect the equivalence class names and PDB IDs, model
# info, and chain info for the representative structures
nrlist_info = parse_nrlist.get_info(nrlist_file)

# collect list of files in mmCIF_files folder
mmCIF_files = []
if os.path.isdir(directory + "/mmCIF_files"):
    mmCIF_files = os.listdir(directory + "/mmCIF_files")
else:
    print("Error: There is no folder named mmCIF_files in the current working directory.")
    sys.exit(1)

# initialize a list that will store dictionaries with information from the representative set file and mmCIF files
eq_class_list = []

# iterate through the equivalence classes identified from the representative set file and gather related information
# from the mmCIF files and store it in dictionaries with the information gathered from the representative set file
for eq_class in nrlist_info:
    # initialize individual variables that store info already collected from the representative set file
    equivalence_class = eq_class[0]
    nrlist_PDB_ID_list = eq_class[1]
    model_list = eq_class[2]
    chain_list = eq_class[3]
    # construct the mmCIF file name
    mmCIF_file_name = nrlist_PDB_ID_list[0].lower() + ".cif"
    if mmCIF_file_name not in mmCIF_files:
        print(f"Error: The mmCIF_files folder does not contain {mmCIF_file_name}.")
        sys.exit(1)
    # read the mmCIF file
    mmCIF = cif.read_file(directory + "/mmCIF_files/" + mmCIF_file_name)
    block = mmCIF.sole_block()
    # collect the PDB ID from the mmCIF file
    mmCIF_PDB_ID = ''
    if block.find_values('_entry.id'):
        pdb = block.find_values('_entry.id').str(0)
        # check if the PDB ID contains the expected number of characters
        if len(pdb) != 4:
            print("Error: A PDB ID from an mmCIF file does not contain the expected number of characters.")
            sys.exit(1)
        mmCIF_PDB_ID = pdb
    else:
        print(f"Warning: {mmCIF_file_name} does not contain _entry.id.")
    # collect the initial release date from the mmCIF file
    initial_release_date = ''
    if block.find_values('_pdbx_audit_revision_history.revision_date'):
        initial_release_date = block.find_values('_pdbx_audit_revision_history.revision_date').str(0)
    else:
        print(f"Warning: {mmCIF_file_name} does not contain _pdbx_audit_revision_history.revision_date.")
    # collect the experimental method from the mmCIF file
    method = ''
    if block.find_values('_exptl.method'):
        method = block.find_values('_exptl.method').str(0)
    else:
        print(f"Warning: {mmCIF_file_name} does not contain _exptl.method.")
    # collect the resolution from the mmCIF file
    resolution = ''
    # resolution for x-ray diffraction structures
    if method == 'X-RAY DIFFRACTION':
        if block.find_values('_reflns.d_resolution_high'):
            resolution = block.find_values('_reflns.d_resolution_high').str(0)
        else:
            print(f"Warning: {mmCIF_file_name} does not contain _reflns.d_resolution_high.")
    # resolution for electron microscopy structures
    elif method == 'ELECTRON MICROSCOPY':
        if block.find_values('_em_3d_reconstruction.resolution'):
            resolution = block.find_values('_em_3d_reconstruction.resolution').str(0)
        else:
            print(f"Warning: {mmCIF_file_name} does not contain _em_3d_reconstruction.resolution.")
    else:
        print(f"Warning: _exptl.method in {mmCIF_file_name} is not 'X-RAY DIFFRACTION' or 'ELECTRON MICROSCOPY'.")
    # collect information related to all the chains in the representative structure
    entity = []
    entity_poly = []
    entity_src_gen = []
    entity_src_nat = []
    pdbx_entity_src_syn = []
    if block.find(['_entity.id', '_entity.src_method', '_entity.pdbx_description']):
        entity = block.find(['_entity.id', '_entity.src_method', '_entity.pdbx_description'])
    else:
        print(f"Error: {mmCIF_file_name} does not contain _entity.id, _entity.src_method, and/or "
              f"_entity.pdbx_description.")
        sys.exit(1)
    if block.find(['_entity_poly.entity_id', '_entity_poly.pdbx_strand_id', '_entity_poly.pdbx_seq_one_letter_code']):
        entity_poly = block.find(['_entity_poly.entity_id', '_entity_poly.pdbx_strand_id',
                                  '_entity_poly.pdbx_seq_one_letter_code', '_entity_poly.pdbx_seq_one_letter_code_can'])
    else:
        print(f"Error: {mmCIF_file_name} does not contain _entity_poly.entity_id, _entity_poly.pdbx_strand_id, "
              f"_entity_poly.pdbx_seq_one_letter_code, and/or _entity_poly.pdbx_seq_one_letter_code_can.")
        sys.exit(1)
    entity_src_exists = False
    if block.find(['_entity_src_gen.entity_id', '_entity_src_gen.pdbx_gene_src_scientific_name']):
        entity_src_gen = block.find(['_entity_src_gen.entity_id', '_entity_src_gen.pdbx_gene_src_scientific_name'])
        entity_src_exists = True
    if block.find(['_entity_src_nat.entity_id', '_entity_src_nat.pdbx_organism_scientific']):
        entity_src_nat = block.find(['_entity_src_nat.entity_id', '_entity_src_nat.pdbx_organism_scientific'])
        entity_src_exists = True
    if block.find(['_pdbx_entity_src_syn.entity_id', '_pdbx_entity_src_syn.organism_scientific']):
        pdbx_entity_src_syn = block.find(['_pdbx_entity_src_syn.entity_id', '_pdbx_entity_src_syn.organism_scientific'])
        entity_src_exists = True
    if entity_src_exists == False:
        print(f"Warning: {mmCIF_file_name} does not contain the designated keywords in any of the categories "
              f"_entity_src_gen, _entity_src_nat, or _pdbx_entity_src_syn.")
    # determine the entity IDs and sequences of the representative chains
    chain_entity_ID_list = []
    sequence_list = []
    sequence_list_can = []
    # iterate through the chains identified in the representative set file
    for eq_class_chain in chain_list:
        find_match = False
        # iterate through the chains in the mmCIF file
        for row in range(len(entity_poly)):
            for mmCIF_chain in entity_poly[row][1].split(','):
                # save the entity ID and sequence of the chain in the mmCIF file that matches the chain in the
                # representative set file
                if eq_class_chain == mmCIF_chain:
                    find_match = True
                    chain_entity_ID_list.append(entity_poly[row][0])
                    sequence_list.append(entity_poly[row][2])
                    sequence_list_can.append(entity_poly[row][3])
        if not find_match:
            print(f"Error: A match between a chain ID in the representative set file and a corresponding chain ID in "
                  f"{mmCIF_file_name} was not found.")
            sys.exit(1)
    # use the chain entity IDs to collect the source and description of each representative chain, then use the chain
    # entity IDs with the source information to collect the organism information for each chain
    chain_src_method_list = []
    chain_description_list = []
    chain_organism_list = []
    for chain_entity_ID in chain_entity_ID_list:
        # iterate through each row of the entity table which contains an entity ID, source information, and a
        # description
        find_match = False
        for row in range(len(entity)):
            # if the entity ID of the chain matches that listed in the row, collect the source and description
            if chain_entity_ID == entity[row][0]:
                find_match = True
                chain_src_method_list.append(entity[row].str(1))
                chain_description_list.append(entity[row].str(2))
        if not find_match:
            chain_src_method_list.append('')
            chain_description_list.append('')
            print(f"Warning: A match between _entity_poly.entity_id and _entity.id in {mmCIF_file_name} was not found.")
        # collect the organism based on the chain entity ID and the source information
        if chain_src_method_list[-1] == 'man':
            # check to ensure entity_src_gen contains information
            if entity_src_gen:
                find_match = False
                for row in range(len(entity_src_gen)):
                    # if the entity ID of the chain matches that listed in the row, collect the organism
                    if chain_entity_ID == entity_src_gen[row][0]:
                        find_match = True
                        chain_organism_list.append(entity_src_gen[row].str(1))
                if not find_match:
                    chain_organism_list.append('')
                    print(f"Warning: A match between _entity.id and _entity_src_gen.entity_id in {mmCIF_file_name} was "
                          f"not found.")
            else:
                chain_organism_list.append('')
                print(f"Warning: {mmCIF_file_name} does not contain _entity_src_gen.entity_id and/or "
                      f"_entity_src_gen.pdbx_gene_src_scientific_name.")
        elif chain_src_method_list[-1] == 'nat':
            # check to ensure entity_src_nat contains information
            if entity_src_nat:
                find_match = False
                for row in range(len(entity_src_nat)):
                    # if the entity ID of the chain matches that listed in the row, collect the organism
                    if chain_entity_ID == entity_src_nat[row][0]:
                        find_match = True
                        chain_organism_list.append(entity_src_nat[row].str(1))
                if not find_match:
                    chain_organism_list.append('')
                    print(f"Warning: A match between _entity.id and _entity_src_nat.entity_id in {mmCIF_file_name} was "
                          f"not found.")
            else:
                chain_organism_list.append('')
                print(f"Warning: {mmCIF_file_name} does not contain _entity_src_nat.entity_id and/or "
                      f"_entity_src_nat.pdbx_organism_scientific.")
        elif chain_src_method_list[-1] == 'syn':
            # check to ensure pdbx_entity_src_syn contains information
            if pdbx_entity_src_syn:
                find_match = False
                for row in range(len(pdbx_entity_src_syn)):
                    # if the entity ID of the chain matches that listed in the row, collect the organism
                    if chain_entity_ID == pdbx_entity_src_syn[row][0]:
                        find_match = True
                        chain_organism_list.append(pdbx_entity_src_syn[row].str(1))
                if not find_match:
                    chain_organism_list.append('')
                    print(f"Warning: A match between _entity.id and _pdbx_entity_src_syn.entity_id in "
                          f"{mmCIF_file_name} was not found.")
            else:
                chain_organism_list.append('')
                print(f"Warning: {mmCIF_file_name} does not contain _pdbx_entity_src_syn.entity_id and/or "
                      f"_pdbx_entity_src_syn.organism_scientific.")
    # create a dictionary with all the gathered information
    eq_class_dict = {'equivalence_class': equivalence_class, 'nrlist_PDB_ID_list': nrlist_PDB_ID_list,
                     'model_list': model_list, 'chain_list': chain_list, 'mmCIF_file_name': mmCIF_file_name,
                     'mmCIF_PDB_ID': mmCIF_PDB_ID, 'initial_release_date': initial_release_date, 'method': method,
                     'resolution': resolution, 'chain_entity_ID_list': chain_entity_ID_list,
                     'sequence_list': sequence_list, 'sequence_list_can': sequence_list_can,
                     'chain_src_method_list': chain_src_method_list, 'chain_description_list': chain_description_list,
                     'chain_organism_list': chain_organism_list}
    # append the dictionary to the list of dictionaries
    eq_class_list.append(eq_class_dict)

# write the information collected to a csv file
keys = ['equivalence_class', 'nrlist_PDB_ID_list', 'model_list', 'chain_list', 'mmCIF_file_name', 'mmCIF_PDB_ID',
        'initial_release_date', 'method', 'resolution', 'chain_entity_ID_list', 'sequence_list', 'sequence_list_can',
        'chain_src_method_list', 'chain_description_list', 'chain_organism_list']
with open("eq_class_info.csv", "w") as csv_file:
    writer = csv.DictWriter(csv_file, fieldnames=keys)
    writer.writeheader()
    writer.writerows(eq_class_list)
