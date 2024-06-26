"""
This script reads information from the equivalence class file created by parse_nrlist.py and specified by the Snakefile.
It works with PyMOL to save an mmCIF file of the original structure to a location specified in the Snakemake
configuration file. Alternative conformations are removed, hydrogens are added, and the modified structure is saved to
an mmCIF file specified by the Snakefile. The script collects data from the modified structure on potential hydrogen
bonds involving atoms of interest and nucleobase atom b-factors in the equivalence class member RNA chains. It writes this data to
three separate csv files specified by the Snakefile. If commit_hash is set to true in the Snakemake configuration file,
the commit hash of the repo will also be written within a commented line to the data files if no uncommitted changes
have been made to the repo.
"""

import sys
import os
import csv
import subprocess
import time
import copy
from pymol import cmd
from pymol import stored
import numpy as np
import pandas as pd
from datetime import datetime
import residue_library
import eval_H_bonding
import remove_alt_conf

# TODO add code to count how many atoms each residue contains pre- and post- hydrogen addition, reject those that do not meet expectations

# Redirect stdout and stderr to log files.
stdout = sys.stdout
stderr = sys.stderr
stdout_file = open(snakemake.log.stdout, mode='w')
stderr_file = open(snakemake.log.stderr, mode='w')
sys.stdout = stdout_file
sys.stderr = stderr_file


# This function is to be run if an error is encountered.
def error(msg_list):

    # Create the output files expected by Snakemake.
    subprocess.run(["touch", snakemake.output.data])

    # If the structure has already been loaded, save the modified structure if specified.
    if cmd.count_atoms('all') > 0 and snakemake.config["save_modified_mmcif"]:
        cmd.save(modified_mmcif_dir + eq_class_mem_id + ".cif")

    # If the structure has already been loaded, remove the original mmCIF if specified.
    if cmd.count_atoms('all') > 0 and snakemake.config["remove_original_mmcif"]:
        subprocess.run(["rm", original_mmcif_dir + f"{pdb_id}.cif".lower()])

    # Print the error message(s).
    for msg in msg_list:
        print(msg)

    # Close files, reset stdout and stderr, and exit.
    stdout_file.close()
    stderr_file.close()
    sys.stdout = stdout
    sys.stderr = stderr
    sys.exit(0)


# Check whether the folder exists for PyMOL to save mmCIF files into from the modified structure.
# If it does not exist, create the folder.
modified_mmcif_dir = snakemake.config["modified_mmcif_dir"]
if not os.path.isdir(modified_mmcif_dir):
    try:
        os.mkdir(modified_mmcif_dir)
    except FileExistsError:
        time.sleep(5.0)
        if not os.path.isdir(modified_mmcif_dir):
            os.mkdir(modified_mmcif_dir)

# If commit_hash is set to true in the Snakemake configuration file, check if any changes have been made to the repo and
# get the hash of the current git commit. If uncommitted changes have been made to anything other than files listed in
# the acceptable_changes variable defined below, print an error message and exit.
repo_changes = []
commit_hash = ""
if snakemake.config["commit_hash"]:
    repo_changes = list(subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=no"],
                                                cwd=os.path.dirname(os.path.realpath(__file__)))
                        .decode('ascii').strip().split("\n"))
    acceptable_changes = ['config/config.yaml', snakemake.config["rep_set_file"], snakemake.config["add_res_file"]]
    for file in repo_changes:
        if file.split(' ')[-1] in acceptable_changes:
            repo_changes.remove(file)
    if len(repo_changes) == 0:
        commit_hash = subprocess.check_output(["git", "rev-parse", "HEAD"],
                                              cwd=os.path.dirname(os.path.realpath(__file__))).decode('ascii').strip()
    else:
        error([f"Error: Uncommitted changes have been made to the repo."])

# Collect the equivalence class name, PDB ID, model info, and chain info from the string describing the equivalence
# class member.
eq_class_build = []
pdb_id_build = []
model_build = []
chain_build = []
index = 0
track = 0
for char in snakemake.wildcards.eq_class_members:
    if index < 3:
        eq_class_build.append(char)
    elif char != '_' and index == 3:
        pdb_id_build.append(char)
    elif char != '_' and index == 4:
        model_build.append(char)
    elif char != '_' and index > 4:
        if index != track:
            chain_build.append([char])
            track = index
        else:
            chain_build[-1].append(char)
    if char == '_':
        index += 1
eq_class = "".join(eq_class_build[:-1])
pdb_id = "".join(pdb_id_build)
model = "".join(model_build)
chain_list = []
for chain in chain_build:
    chain_list.append("".join(chain))

# Create an identifier for the equivalence class member.
eq_class_mem_id = f'{eq_class}_{pdb_id}_{model}{"".join(["_" + chain for chain in chain_list])}'

# Check whether the folder exists for PyMOL to save mmCIF files into that are fetched from the PDB.
# If it does not exist, create the folder.
original_mmcif_dir = snakemake.config["original_mmcif_dir"]
if not os.path.isdir(original_mmcif_dir):
    try:
        os.mkdir(original_mmcif_dir)
    except FileExistsError:
        time.sleep(5.0)
        if not os.path.isdir(original_mmcif_dir):
            os.mkdir(original_mmcif_dir)

# Change fetch_path to original_mmcif_dir so that PyMOL drops the mmCIF files into this folder when running fetch.
cmd.set('fetch_path', cmd.exp_path(original_mmcif_dir))

# Construct two strings that can be used with PyMOL to select donor and rotatable donor atoms.
donor_string = ''
rotatable_donor_string = ''
for residue in residue_library.RESIDUE_LIBRARY:
    for donor in residue['don']:
        donor_string += f'(resn {residue["res"]} and name {donor[0]}) '
        if donor[2]:
            rotatable_donor_string += f'(resn {residue["res"]} and name {donor[0]}) '
donor_string = donor_string[:-1]
rotatable_donor_string = rotatable_donor_string[:-1]

# Construct a string that can be used with PyMOL to select all possible donor atoms of particular interest.
donors_of_interest_str = ''
for atom in snakemake.config["donors_of_interest"]:
    donors_of_interest_str += f'(resn {atom.split(".")[0]} and name {atom.split(".")[1]}) '
donors_of_interest_str = donors_of_interest_str[:-1]

# Prepare a string that can be used with PyMOL to identify the chains of all the equivalence class member RNAs.
mem_rna_chains = " ".join(["chain " + chain for chain in chain_list])

# Retrieve the structure that contains the equivalence class member RNA chains.
cmd.fetch(pdb_id)

# If the structure contains more than one model (or "state" in PyMOL), create a new PyMOL object of just that model.
if cmd.count_states(pdb_id) > 1:
    cmd.create(f'{pdb_id}_state_{model}', selection=pdb_id,
               source_state=model,
               target_state=1)
    cmd.delete(pdb_id)

# Remove atoms representing alternative conformations.
remove_status = remove_alt_conf.remove(eq_class_mem_id)
successful_completion = remove_status[0]
if not successful_completion:
    error(remove_status[1])

# Remove any hydrogens that loaded with the structure.
cmd.remove('elem H')

# Add hydrogens to non-rotatable donors that are a part of or near the equivalence class member RNA.
cmd.h_add(f'(({donor_string}) and not ({rotatable_donor_string})) within {snakemake.config["search_dist"]} of '
          f'({mem_rna_chains})')

# Randomly sample 100 atom indices in the structure and record the info of the associated atoms for later comparison.
num_atoms = cmd.count_atoms('all')
indices = np.random.default_rng().integers(low=1, high=num_atoms+1, size=100)
stored.check_one = []
for index in indices:
    cmd.iterate(f'index {index}', 'stored.check_one.append([index, name, resn, resi, chain])')

# Store a list of donors of interest from the equivalence class member RNA chains.
stored.donor_list = []
cmd.iterate(f'({mem_rna_chains}) and ({donors_of_interest_str})',
            'stored.donor_list.append([index, name, resn, resi, chain, segi])')

# Store the number of heavy atoms belonging to organic molecules or polymers near the donor of interest within a
# dictionary. Additionally, store the average b-factor of the heavy atoms that make up the nucleobase containing the
# donors of interest. Values of one for the DOI key indicate that these are donors of interest.
info_dict = {"don_index": [], "don_name": [], "don_resn": [], "don_resi": [], "don_chain": [], "don_segi": [],
             "count_1": [], "count_2": [], "b-factor": [], "DOI": []}
for donor in stored.donor_list:
    # Count the number of heavy atoms belonging to organic molecules or polymers near the donor. Exclude the nucleobase
    # of the donor.
    count_1 = cmd.count_atoms(f'((organic or polymer) within {snakemake.config["count_dist_1"]} of index {donor[0]}) '
                              f'and not ((sidechain and byres index {donor[0]}) or elem H)')
    count_2 = cmd.count_atoms(f'((organic or polymer) within {snakemake.config["count_dist_2"]} of index {donor[0]}) '
                              f'and not ((sidechain and byres index {donor[0]}) or elem H)')
    # Collect the b-factors for side chain atoms.
    stored.b_factors = []
    cmd.iterate(f"resn {donor[2]} and resi {donor[3]} and chain {donor[4]} and sidechain and not elem H",
                "stored.b_factors.append(b)")
    # Calculate the average of the b-factors.
    b_factor_avg = sum(stored.b_factors)/len(stored.b_factors)
    # Add the donor, heavy atom counts, and average b-factor to the count dictionary.
    info_row = donor + [count_1, count_2, b_factor_avg, 1]
    for key, value in zip(info_dict, info_row):
        info_dict[key].append(value)

# Create a dataframe based on the master dictionary.
info_df = pd.DataFrame(info_dict)

# Store a list of acceptors near the donors of interest within a dictionary of nucleobases.
atom_pair_dict = {"don_index": [], "don_name": [], "don_resn": [], "don_resi": [], "don_chain": [], "don_segi": [],
                  "acc_index": [], "acc_name": [], "acc_resn": [], "acc_resi": [], "acc_chain": [], "acc_segi": []}
for donor in stored.donor_list:
    # Find the acceptors near the donor. Exclude the nucleobase of the donor.
    stored.acceptors = []
    cmd.iterate(f'(acceptors within {snakemake.config["search_dist"]} of index {donor[0]}) and (organic or polymer) '
                f'and not (sidechain and byres index {donor[0]})',
                'stored.acceptors.append([index, name, resn, resi, chain, segi])')
    # Add the acceptors to the nucleobase dictionary.
    for acceptor in stored.acceptors:
        atom_pair_row = donor + acceptor
        for key, value in zip(atom_pair_dict, atom_pair_row):
            atom_pair_dict[key].append(value)

# Create a dataframe based on the nucleobase dictionary.
atom_pair_df = pd.DataFrame(atom_pair_dict)

# Acquire the H-bonding geometry measurements for all acceptors near each donor of interest. Values of one for the DOI
# key indicate that the donor atoms are donors of interest.
don_h_bonds_dict = {"don_index": [], "don_name": [], "don_resn": [], "don_resi": [], "don_chain": [], "don_segi": [],
                    "acc_index": [], "acc_name": [], "acc_resn": [], "acc_resi": [], "acc_chain": [], "acc_segi": [],
                    "don_acc_distance": [], "h_acc_distance": [], "h_angle": [], "h_dihedral": [], "h_name": [],
                    "DOI": []}
for atom_pair in atom_pair_df.itertuples():
    # Store the donor and acceptor atom values for the dataframe row.
    don_list = [atom_pair.don_index, atom_pair.don_name, atom_pair.don_resn, atom_pair.don_resi, atom_pair.don_chain]
    acc_list = [atom_pair.acc_index, atom_pair.acc_name, atom_pair.acc_resn, atom_pair.acc_resi, atom_pair.acc_chain]
    # Retrieve the H-bond measurements for the atom pair.
    h_bond_list = eval_H_bonding.evaluate(don_list, acc_list, eq_class_mem_id, residue_library.RESIDUE_LIBRARY)
    # If the H-bond evaluation is not successful, print the error message(s) and exit.
    successful_completion = h_bond_list[0]
    if not successful_completion:
        error(h_bond_list[1])
    # Add the H-bond measurements to the dictionary.
    else:
        for h_bond in h_bond_list[1]:
            row = don_list + [atom_pair.don_segi] + acc_list + [atom_pair.acc_segi] + h_bond + [1]
            for key, value in zip(don_h_bonds_dict, row):
                don_h_bonds_dict[key].append(value)

# Create a dataframe based on the dictionary.
don_h_bonds_df = pd.DataFrame(don_h_bonds_dict).astype({"h_dihedral": "object"})

# Store a list of G(N1) donors (or its deoxy counterpart) from within or near the equivalence class member RNA chains.
stored.g_n1_list = []
cmd.iterate(f'({mem_rna_chains} expand {snakemake.config["search_dist"]}) and ((resn G or resn DG) and name N1)',
            'stored.g_n1_list.append([index, name, resn, resi, chain, segi])')

# Store a list of atoms near the G(N1) donors within a G(N1) dictionary.
g_n1_dict = {"don_index": [], "don_name": [], "don_resn": [], "don_resi": [], "don_chain": [], "don_segi": [],
             "acc_index": [], "acc_name": [], "acc_resn": [], "acc_resi": [], "acc_chain": [], "acc_segi": []}
for donor in stored.g_n1_list:
    # Find the atoms near the donor.
    stored.nearby_atoms = []
    cmd.iterate(f'index {donor[0]} around {snakemake.config["search_dist"]}',
                'stored.nearby_atoms.append([index, name, resn, resi, chain, segi])')
    # If the nearby atom is a C(N3), U(O2), or their deoxy counterparts, add them to the G(N1) dictionary along with the
    # info for the G(N1) donor.
    for atom in stored.nearby_atoms:
        if (atom[2], atom[1]) in [("C", "N3"), ("U", "O2"), ("DC", "N3"), ("DT", "O2")]:
            g_n1_row = donor + atom
            for key, value in zip(g_n1_dict, g_n1_row):
                g_n1_dict[key].append(value)

# Create a dataframe based on the G(N1) dictionary.
g_n1_df = pd.DataFrame(g_n1_dict)

# Acquire the H-bonding geometry measurements for the acceptors near each G(N1). Values of zero for the DOI key indicate
# that the donor atoms are not donors of interest.
g_n1_h_bonds_dict = {"don_index": [], "don_name": [], "don_resn": [], "don_resi": [], "don_chain": [], "don_segi": [],
                     "acc_index": [], "acc_name": [], "acc_resn": [], "acc_resi": [], "acc_chain": [], "acc_segi": [],
                     "don_acc_distance": [], "h_acc_distance": [], "h_angle": [], "h_dihedral": [], "h_name": [],
                     "DOI": []}
for atom_pair in g_n1_df.itertuples():
    # Store the donor and acceptor atom values for the dataframe row.
    don_list = [atom_pair.don_index, atom_pair.don_name, atom_pair.don_resn, atom_pair.don_resi, atom_pair.don_chain]
    acc_list = [atom_pair.acc_index, atom_pair.acc_name, atom_pair.acc_resn, atom_pair.acc_resi, atom_pair.acc_chain]
    # Retrieve the H-bond measurements for the atom pair.
    h_bond_list = eval_H_bonding.evaluate(don_list, acc_list, eq_class_mem_id, residue_library.RESIDUE_LIBRARY)
    # If the H-bond evaluation is not successful, print the error message(s) and exit.
    successful_completion = h_bond_list[0]
    if not successful_completion:
        error(h_bond_list[1])
    # Add the H-bond measurements to the dictionary.
    else:
        for h_bond in h_bond_list[1]:
            row = don_list + [atom_pair.don_segi] + acc_list + [atom_pair.acc_segi] + h_bond + [0]
            for key, value in zip(g_n1_h_bonds_dict, row):
                g_n1_h_bonds_dict[key].append(value)

# Create a dataframe based on the dictionary.
g_n1_h_bonds_df = pd.DataFrame(g_n1_h_bonds_dict).astype({"h_dihedral": "object"})

# Store a list of U(N3) donors (or its deoxy counterpart) from within or near the equivalence class member RNA chains.
stored.u_n3_list = []
cmd.iterate(f'({mem_rna_chains} expand {snakemake.config["search_dist"]}) and ((resn U or resn DT) and name N3)',
            'stored.u_n3_list.append([index, name, resn, resi, chain, segi])')

# Store a list of atoms near the U(N3) donors within a U(N3) dictionary.
u_n3_dict = {"don_index": [], "don_name": [], "don_resn": [], "don_resi": [], "don_chain": [], "don_segi": [],
             "acc_index": [], "acc_name": [], "acc_resn": [], "acc_resi": [], "acc_chain": [], "acc_segi": []}
for donor in stored.u_n3_list:
    # Find the atoms near the donor.
    stored.nearby_atoms = []
    cmd.iterate(f'index {donor[0]} around {snakemake.config["search_dist"]}',
                'stored.nearby_atoms.append([index, name, resn, resi, chain, segi])')
    # If the nearby atom is a A(N1), G(O6), or their deoxy counterparts, add them to the U(N3) dictionary along with the
    # info for the U(N3) donor.
    for atom in stored.nearby_atoms:
        if (atom[2], atom[1]) in [("A", "N1"), ("G", "O6"), ("DA", "N1"), ("DG", "O6")]:
            u_n3_row = donor + atom
            for key, value in zip(u_n3_dict, u_n3_row):
                u_n3_dict[key].append(value)

# Create a dataframe based on the U(N3) dictionary.
u_n3_df = pd.DataFrame(u_n3_dict)

# Acquire the H-bonding geometry measurements for the acceptors near each U(N3). Values of zero for the DOI key indicate
# that the donor atoms are not donors of interest.
u_n3_h_bonds_dict = {"don_index": [], "don_name": [], "don_resn": [], "don_resi": [], "don_chain": [], "don_segi": [],
                     "acc_index": [], "acc_name": [], "acc_resn": [], "acc_resi": [], "acc_chain": [], "acc_segi": [],
                     "don_acc_distance": [], "h_acc_distance": [], "h_angle": [], "h_dihedral": [], "h_name": [],
                     "DOI": []}
for atom_pair in u_n3_df.itertuples():
    # Store the donor and acceptor atom values for the dataframe row.
    don_list = [atom_pair.don_index, atom_pair.don_name, atom_pair.don_resn, atom_pair.don_resi, atom_pair.don_chain]
    acc_list = [atom_pair.acc_index, atom_pair.acc_name, atom_pair.acc_resn, atom_pair.acc_resi, atom_pair.acc_chain]
    # Retrieve the H-bond measurements for the atom pair.
    h_bond_list = eval_H_bonding.evaluate(don_list, acc_list, eq_class_mem_id, residue_library.RESIDUE_LIBRARY)
    # If the H-bond evaluation is not successful, print the error message(s) and exit.
    successful_completion = h_bond_list[0]
    if not successful_completion:
        error(h_bond_list[1])
    # Add the H-bond measurements to the dictionary.
    else:
        for h_bond in h_bond_list[1]:
            row = don_list + [atom_pair.don_segi] + acc_list + [atom_pair.acc_segi] + h_bond + [0]
            for key, value in zip(u_n3_h_bonds_dict, row):
                u_n3_h_bonds_dict[key].append(value)

# Create a dataframe based on the dictionary.
u_n3_h_bonds_df = pd.DataFrame(u_n3_h_bonds_dict).astype({"h_dihedral": "object"})

# Record info on the atoms associated with the previously determined 100 atom indices and check whether there are any
# changes.
stored.check_two = []
for index in indices:
    cmd.iterate(f'index {index}', 'stored.check_two.append([index, name, resn, resi, chain])')
if not stored.check_one == stored.check_two:
    error([f"Error: The indices in equivalence class member {snakemake.wildcards.eq_class_members} changed."])

# Prepare a master dataframe containing heavy atom count, b-factor, H-bonding data, and other relevant information.
master_df = info_df.merge(don_h_bonds_df, how='outer').merge(g_n1_h_bonds_df, how='outer').merge(u_n3_h_bonds_df,
                                                                                                 how='outer')
master_df.loc[:, ['model', 'PDB', 'eq_class_member']] = [model, pdb_id, eq_class_mem_id]

# Write a csv containing the data stored in the master dataframe.
with open(snakemake.output.data, "w") as csv_file:
    writer = csv.writer(csv_file)
    if commit_hash:
        writer.writerow([f"# dual-H-bonding-nucleobases repo git commit hash: {commit_hash}"])
    writer.writerow([f"# representative set file: {snakemake.config['rep_set_file']}"])
    writer.writerow([f"# file created on: {datetime.now().strftime('%y-%m-%d %H:%M:%S.%f')}"])
master_df.to_csv(snakemake.output.data, index=False, mode='a', na_rep='NaN')

# Save the modified structure if specified.
if snakemake.config["save_modified_mmcif"]:
    cmd.save(modified_mmcif_dir + eq_class_mem_id + ".cif")

# Remove the original mmCIF if specified.
if snakemake.config["remove_original_mmcif"]:
    subprocess.run(["rm", original_mmcif_dir + f"{pdb_id}.cif".lower()])

# Close files and reset stdout and stderr.
stdout_file.close()
stderr_file.close()
sys.stdout = stdout
sys.stderr = stderr
