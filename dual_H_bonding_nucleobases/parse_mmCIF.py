"""
This module reads the representative set file and iterates through the mmCIF files to collect information related to
each equivalence class that may be useful for data analysis.
"""

import sys
import os
import csv
from gemmi import cif

directory = os.getcwd()

# initialize a dictionary that will store information from the representative set file and mmCIF files
eq_class_info = {'equivalence_class': [], 'nrlist_PDB_ID': [], 'model': [], 'chain': [], 'mmCIF_file_name': [],
                 'mmCIF_PDB_ID': [], 'initial_release_date': [], 'method': [], 'resolution': [], 'chain_entity_ID': [],
                 'chain_organism': [], 'chain_src_method': [], 'chain_description': []}

# identify the representative set file
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
        nrlist_info.append([line[0], PDB_ID_list, model_list, chain_list])

# collect list of files in mmCIF_files folder
mmCIF_files = []
if os.path.isdir(directory + "/mmCIF_files"):
    mmCIF_files = os.listdir(directory + "/mmCIF_files")
else:
    print("Error: There is no folder named mmCIF_files in the current working directory.")
    sys.exit(1)

# collect information from the mmCIF files
for i in range(len(nrlist_info)):
    eq_class_info['equivalence_class']
    eq_class_info['nrlist_PDB_ID']
    eq_class_info['model']
    eq_class_info['chain']
    mmCIF_name = eq_class_info['nrlist_PDB_ID'][i][0].lower() + ".cif"
    if mmCIF_name not in mmCIF_files:
        print(f"Error: The mmCIF_files folder does not contain {mmCIF_name}.")
        sys.exit(1)
    eq_class_info['mmCIF_file_name'].append(mmCIF_name)
    mmCIF = cif.read_file(directory + "/mmCIF_files/" + mmCIF_name)
    block = mmCIF.sole_block()
    # collect the PDB ID from the mmCIF file
    if block.find_values('_entry.id'):
        pdb = block.find_values('_entry.id').str(0)
        # check if the PDB ID contains the expected number of characters
        if len(pdb) != 4:
            print("Error: A PDB ID from an mmCIF file does not contain the expected number of characters.")
            sys.exit(1)
        eq_class_info['mmCIF_PDB_ID'].append(pdb)
    else:
        eq_class_info['mmCIF_PDB_ID'].append("")
        print(f"Warning: {mmCIF_name} does not contain _entry.id.")
    # collect the initial release date from the mmCIF file
    if block.find_values('_pdbx_audit_revision_history.revision_date'):
        release = block.find_values('_pdbx_audit_revision_history.revision_date').str(0)
        eq_class_info['initial_release_date'].append(release)
    else:
        eq_class_info['initial_release_date'].append("")
        print(f"Warning: {mmCIF_name} does not contain _pdbx_audit_revision_history.revision_date.")
    # collect the experimental method from the mmCIF file
    if block.find_values('_exptl.method'):
        method = block.find_values('_exptl.method').str(0)
        eq_class_info['method'].append(method)
    else:
        eq_class_info['method'].append("")
        print(f"Warning: {mmCIF_name} does not contain _exptl.method.")
    # collect the resolution from the mmCIF file
    # resolution for x-ray diffraction structures
    if method == 'X-RAY DIFFRACTION':
        if block.find_values('_reflns.d_resolution_high'):
            resolution = block.find_values('_reflns.d_resolution_high').str(0)
            eq_class_info['resolution'].append(resolution)
        else:
            eq_class_info['resolution'].append("")
            print(f"Warning: {mmCIF_name} does not contain _reflns.d_resolution_high.")
    # resolution for electron microscopy structures
    elif method == 'ELECTRON MICROSCOPY':
        if block.find_values('_em_3d_reconstruction.resolution'):
            resolution = block.find_values('_em_3d_reconstruction.resolution').str(0)
            eq_class_info['resolution'].append(resolution)
        else:
            eq_class_info['resolution'].append("")
            print(f"Warning: {mmCIF_name} does not contain _em_3d_reconstruction.resolution.")
    else:
        eq_class_info['resolution'].append("")
    # collect information related to all the chains in the representative structure of each equivalence class from the
    # mmCIF file
    entity = []
    entity_poly = []
    entity_src_gen = []
    entity_src_nat = []
    pdbx_entity_src_syn = []
    if block.find(['_entity.id', '_entity.src_method', '_entity.pdbx_description']):
        entity = block.find(['_entity.id', '_entity.src_method', '_entity.pdbx_description'])
    else:
        print(f"Error: {mmCIF_name} does not contain _entity.id, _entity.src_method, and/or "
              f"_entity.pdbx_description.")
        sys.exit(1)
    if block.find(['_entity_poly.entity_id', '_entity_poly.pdbx_strand_id']):
        entity_poly = block.find(['_entity_poly.entity_id', '_entity_poly.pdbx_strand_id'])
    else:
        print(f"Error: {mmCIF_name} does not contain _entity_poly.entity_id and/or _entity_poly.pdbx_strand_id.")
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
        print(f"Warning: {mmCIF_name} does not contain the designated keywords in any of the categories "
              f"_entity_src_gen, _entity_src_nat, or _pdbx_entity_src_syn.")
    # determine the entity IDs of the representative chains
    chain_entity_id_list = []
    # iterate through the chains identified in the representative set file
    for eq_class_chain in eq_class_info['chain'][i]:
        # iterate through the chains in the mmCIF file
        for row in range(len(entity_poly)):
            for mmCIF_chain in entity_poly[row][1].split(','):
                # save the entity ID of the chain in the mmCIF file that matches the chain in the representative set
                # file
                if eq_class_chain == mmCIF_chain:
                    chain_entity_id_list.append(entity_poly[row][0])
    eq_class_info['chain_entity_ID'].append(chain_entity_id_list)
    # use the chain entity IDs to collect the source and description of each representative chain, then use the chain
    # entity IDs with the source information to collect the organism information for each chain
    chain_src_method_list = []
    chain_description_list = []
    chain_organism_list = []
    for chain_entity_id in chain_entity_id_list:
        # iterate through each row of the entity table which contains an entity ID, source information, and a
        # description
        for row in range(len(entity)):
            # if the entity ID of the chain matches that listed in the row, collect the source and description
            if chain_entity_id == entity[row][0]:
                chain_src_method_list.append(entity[row].str(1))
                chain_description_list.append(entity[row].str(2))
        # collect the organism based on the chain entity ID and the source information
        if chain_src_method_list[-1] == 'man':
            # check to ensure entity_src_gen contains information
            if entity_src_gen:
                for row in range(len(entity_src_gen)):
                    # if the entity ID of the chain matches that listed in the row, collect the organism
                    if chain_entity_id == entity_src_gen[row][0]:
                        chain_organism_list.append(entity_src_gen[row].str(1))
            else:
                chain_organism_list.append("")
                print(f"Warning: {mmCIF_name} does not contain _entity_src_gen.entity_id and/or "
                      f"_entity_src_gen.pdbx_gene_src_scientific_name.")
        elif chain_src_method_list[-1] == 'nat':
            # check to ensure entity_src_nat contains information
            if entity_src_nat:
                for row in range(len(entity_src_nat)):
                    # if the entity ID of the chain matches that listed in the row, collect the organism
                    if chain_entity_id == entity_src_nat[row][0]:
                        chain_organism_list.append(entity_src_nat[row].str(1))
            else:
                chain_organism_list.append("")
                print(f"Warning: {mmCIF_name} does not contain _entity_src_nat.entity_id and/or "
                      f"_entity_src_nat.pdbx_organism_scientific.")
        elif chain_src_method_list[-1] == 'syn':
            # check to ensure pdbx_entity_src_syn contains information
            if pdbx_entity_src_syn:
                for row in range(len(pdbx_entity_src_syn)):
                    # if the entity ID of the chain matches that listed in the row, collect the organism
                    if chain_entity_id == pdbx_entity_src_syn[row][0]:
                        chain_organism_list.append(pdbx_entity_src_syn[row].str(1))
            else:
                chain_organism_list.append("")
                print(f"Warning: {mmCIF_name} does not contain _pdbx_entity_src_syn.entity_id and/or "
                      f"_pdbx_entity_src_syn.organism_scientific.")
    eq_class_info['chain_src_method'].append(chain_src_method_list)
    eq_class_info['chain_description'].append(chain_description_list)
    eq_class_info['chain_organism'].append(chain_organism_list)

# check to ensure that all the values in the eq_class_info dictionary are all equal length
if not (len(eq_class_info['equivalence_class']) == len(eq_class_info['nrlist_PDB_ID']) == len(eq_class_info['model']) ==
        len(eq_class_info['chain']) == len(eq_class_info['mmCIF_file_name']) == len(eq_class_info['mmCIF_PDB_ID']) ==
        len(eq_class_info['initial_release_date']) == len(eq_class_info['method']) ==
        len(eq_class_info['resolution']) == len(eq_class_info['chain_entity_ID']) ==
        len(eq_class_info['chain_organism']) == len(eq_class_info['chain_src_method']) ==
        len(eq_class_info['chain_description'])):
    print("Error: The values in the eq_class_info dictionary should all be equal length.")
    sys.exit(1)

# with open("eq_class_info.csv", mode='w') as file:
#     for i in range(len(eq_class_info['equivalence_class'])):

print(eq_class_info.values())


#         fields = ['equivalence_class', 'nrlist_PDB_ID', 'model', 'chain', 'mmCIF_file_name', 'mmCIF_PDB_ID',
#                   'initial_release_date', 'method', 'resolution', 'chain_entity_ID', 'chain_organism',
#                   'chain_src_method', 'chain_description']