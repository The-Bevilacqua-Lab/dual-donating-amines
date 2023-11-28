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
                 'mmCIF_PDB_ID': [], 'initial_release_date': [], 'method': [], 'resolution': [], 'organism': []}

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
        eq_class_info['equivalence_class'].append(line[0])
        eq_class_info['nrlist_PDB_ID'].append(PDB_ID_list)
        eq_class_info['model'].append(model_list)
        eq_class_info['chain'].append(chain_list)

# collect list of files in mmCIF_files folder
mmCIF_files = []
if os.path.isdir(directory + "/mmCIF_files"):
    mmCIF_files = os.listdir(directory + "/mmCIF_files")
else:
    print("Error: There is no folder named mmCIF_files in the current working directory.")
    sys.exit(1)

# collect information from the mmCIF files
for i in range(len(eq_class_info['equivalence_class'])):
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
    # collect the initial release date from the mmCIF file
    if block.find_values('_pdbx_audit_revision_history.revision_date'):
        release = block.find_values('_pdbx_audit_revision_history.revision_date').str(0)
        eq_class_info['initial_release_date'].append(release)
    else:
        eq_class_info['initial_release_date'].append("")
    # collect the experimental method from the mmCIF file
    if block.find_values('_exptl.method'):
        method = block.find_values('_exptl.method').str(0)
        eq_class_info['method'].append(method)
    else:
        eq_class_info['method'].append("")
    # collect the resolution from the mmCIF file
    # resolution for x-ray diffraction structures
    if method == 'X-RAY DIFFRACTION':
        if block.find_values('_reflns.d_resolution_high'):
            resolution = block.find_values('_reflns.d_resolution_high').str(0)
            eq_class_info['resolution'].append(resolution)
        else:
            eq_class_info['resolution'].append("")
    # resolution for electron microscopy structures
    elif method == 'ELECTRON MICROSCOPY':
        if block.find_values('_em_3d_reconstruction.resolution'):
            resolution = block.find_values('_em_3d_reconstruction.resolution').str(0)
            eq_class_info['resolution'].append(resolution)
        else:
            eq_class_info['resolution'].append("")
    else:
        eq_class_info['resolution'].append("")
    # collect information related to the individual representative chains of each equivalence class from the mmCIF file
    if block.find(['_entity.id', '_entity.src_method', '_entity.pdbx_description']):
        entity = block.find(['_entity.id', '_entity.src_method', '_entity.pdbx_description'])
    if block.find(['_entity_poly.entity_id', '_entity_poly.pdbx_strand_id']):
        entity_poly = block.find(['_entity_poly.entity_id', '_entity_poly.pdbx_strand_id'])
    if block.find(['_entity_src_gen.entity_id', '_entity_src_gen.pdbx_gene_src_scientific_name']):
        entity_src_gen
    if block.find(['_entity_src_nat.entity_id', '_entity_src_nat.pdbx_organism_scientific']):
        entity_src_nat
    if block.find(['_pdbx_entity_src_syn.entity_id', '_pdbx_entity_src_syn.organism_scientific']):
        pdbx_entity_src_syn



#     # collect the organism from the mmCIF file
#     # organism for x-ray diffraction structures
#     if method == 'X-RAY DIFFRACTION':
#         if block.find_values('_pdbx_entity_src_syn.organism_scientific'):
#             organism = block.find_values('_pdbx_entity_src_syn.organism_scientific').str(0)
#             eq_class_info['organism'].append(organism)
#         else:
#             eq_class_info['organism'].append("")
#     # organism for electron microscopy structures
#     elif method == 'ELECTRON MICROSCOPY':
#         if block.find_values('_em_entity_assembly_naturalsource.organism'):
#             organism = block.find_values('_em_entity_assembly_naturalsource.organism').str(0)
#             eq_class_info['organism'].append(organism)
#         else:
#             eq_class_info['organism'].append("")
#     else:
#         eq_class_info['organism'].append("")
#
# print(eq_class_info['organism'])
# for f in mmCIF_files:
#     if f[-4:] == ".cif":
#         print(f)
