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
        print(f"Error: Uncommitted changes have been made to the repo.")
        sys.exit(1)

# Iterate through the lines of the equivalence class file and collect the equivalence class name, PDB ID, model info,
# and chain info for the equivalence class member.
with open(snakemake.input[0], mode='r') as read_file:
    eq_class_mem = []
    for line in csv.reader(read_file):
        if line[0][0] != "#":
            eq_class_mem.append(line)

# Create an identifier for the equivalence class member.
eq_class_mem_id = (f'{eq_class_mem[0][0]}_{eq_class_mem[1][0]}_{eq_class_mem[2][0]}_'
                   f'{"".join([chain + "_" for chain in eq_class_mem[3]])}')[:-1]

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

# Define the donor and acceptor atoms of the RNA nucleobases.
nuc_donors = ["A.N6", "C.N4", "G.N1", "G.N2", "U.N3"]
nuc_acceptors = ["A.N1", "A.N3", "A.N7", "C.O2", "C.N3", "G.N3", "G.O6", "G.N7", "U.O2", "U.O4"]

# Identify the chains of all the equivalence class member RNAs.
rna_list = []
for idx in range(len(eq_class_mem[2])):
    rna_list.append(f"chain {eq_class_mem[3][idx]}")
mem_rna_chains = " ".join(rna_list)

# Retrieve the structure that contains the equivalence class member RNA chains.
cmd.fetch(eq_class_mem[1][0])

# If the structure contains more than one model (or "state" in PyMOL), create a new PyMOL object of just that model.
if cmd.count_states(eq_class_mem[1][0]) > 1:
    cmd.create(f'{eq_class_mem[1][0]}_state_{eq_class_mem[2][0]}', selection=eq_class_mem[1][0],
               source_state=eq_class_mem[2][0],
               target_state=1)
    cmd.delete(eq_class_mem[1][0])

# Remove atoms representing alternative conformations.
remove_status = remove_alt_conf.remove(eq_class_mem_id)
successful_completion = remove_status[0]
if not successful_completion:
    for note in remove_status[1]:
        print(note)
    sys.exit(1)

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
    cmd.iterate(f'index {index}', 'stored.check_one.append((index, name, resn, resi, chain))')

# Store a list of donors of interest from the equivalence class member RNA chains.
stored.donor_list = []
cmd.iterate(f'({mem_rna_chains}) and ({donors_of_interest_str})',
            'stored.donor_list.append((index, name, resn, resi, chain))')

# Store a list of distances for heavy atoms near the donors of interest within a dictionary.
distance_dict = {"index": [], "name": [], "resn": [], "resi": [], "chain": [], "distance": []}
for donor in stored.donor_list:
    # Get the distances of heavy atoms belonging to organic molecules or polymers near the donor. Exclude the nucleobase
    # of the donor.
    stored.distances = []
    stored.donor_index = donor[0]
    cmd.iterate(f'((organic or polymer) within {snakemake.config["density_dist"]} of index {stored.donor_index}) and '
                f'not ((sidechain and byres index {stored.donor_index}) or elem H)',
                'stored.distances.append(cmd.get_distance(f"index {stored.donor_index}", f"index {index}"))')
    # Add the donor and corresponding distances to the distance dictionary.
    for distance in stored.distances:
        distance_dict["index"].append(donor[0])
        distance_dict["name"].append(donor[1])
        distance_dict["resn"].append(donor[2])
        distance_dict["resi"].append(donor[3])
        distance_dict["chain"].append(donor[4])
        distance_dict["distance"].append(distance)

# Create a dataframe based on the distance dictionary.
distance_df = pd.DataFrame(distance_dict)

# Store a list of acceptors near the donors of interest within a dictionary of nucleobases. Also include two keys which
# have values containing both residue and atom names for donor and acceptor atoms.
nucleobase_dict = {"don_index": [], "don_name": [], "don_resn": [], "don_resi": [], "don_chain": [],
                   "acc_index": [], "acc_name": [], "acc_resn": [], "acc_resi": [], "acc_chain": []}
for donor in stored.donor_list:
    # Find the acceptors near the donor. Exclude the nucleobase of the donor.
    stored.acceptors = []
    cmd.iterate(f'(acceptors within {snakemake.config["search_dist"]} of index {donor[0]}) and (organic or polymer) '
                f'and not (sidechain and byres index {donor[0]})',
                'stored.acceptors.append((index, name, resn, resi, chain))')
    for acceptor in stored.acceptors:
        # Add the donor to the nucleobase dictionary.
        nucleobase_dict["don_index"].append(donor[0])
        nucleobase_dict["don_name"].append(donor[1])
        nucleobase_dict["don_resn"].append(donor[2])
        nucleobase_dict["don_resi"].append(donor[3])
        nucleobase_dict["don_chain"].append(donor[4])
        # Add the acceptors to the nucleobase dictionary.
        nucleobase_dict["acc_index"].append(acceptor[0])
        nucleobase_dict["acc_name"].append(acceptor[1])
        nucleobase_dict["acc_resn"].append(acceptor[2])
        nucleobase_dict["acc_resi"].append(acceptor[3])
        nucleobase_dict["acc_chain"].append(acceptor[4])

# Create a dataframe based on the nucleobase dictionary.
nucleobase_df = pd.DataFrame(nucleobase_dict)

# Acquire the H-bonding geometry measurements for all acceptors near each donor of interest.
don_h_bonds_dict = {"don_index": [], "don_name": [], "don_resn": [], "don_resi": [], "don_chain": [], "acc_index": [],
                    "acc_name": [], "acc_resn": [], "acc_resi": [], "acc_chain": [], "don_acc_distance": [],
                    "h_acc_distance": [], "h_angle": [], "h_dihedral": [], "h_name": []}
for atom_pair in nucleobase_df.itertuples():
    # Store the donor and acceptor atom values for the dataframe row.
    don_list = [atom_pair.don_index, atom_pair.don_name, atom_pair.don_resn, atom_pair.don_resi, atom_pair.don_chain]
    acc_list = [atom_pair.acc_index, atom_pair.acc_name, atom_pair.acc_resn, atom_pair.acc_resi, atom_pair.acc_chain]
    # Retrieve the H-bond measurements for the atom pair.
    h_bond_list = eval_H_bonding.evaluate(don_list, acc_list, eq_class_mem_id, residue_library.RESIDUE_LIBRARY)
    # If the H-bond evaluation is not successful, print the error message(s) and exit.
    successful_completion = h_bond_list[0]
    if not successful_completion:
        for note in h_bond_list[1]:
            print(note)
        sys.exit(1)
    # Add the H-bond measurements to the dictionary.
    else:
        for h_bond in h_bond_list[1]:
            row = don_list + acc_list + h_bond
            idx = 0
            for key in don_h_bonds_dict:
                don_h_bonds_dict[key].append(row[idx])
                idx += 1
# Create a dataframe based on the dictionary.
don_h_bonds_df = pd.DataFrame(don_h_bonds_dict)

# Store a list of G(N1) donors (or its deoxy counterpart) from within or near the equivalence class member RNA chains.
stored.g_n1_list = []
cmd.iterate(f'({mem_rna_chains} expand {snakemake.config["search_dist"]}) and ((resn G or resn DG) and name N1)',
            'stored.g_n1_list.append((index, name, resn, resi, chain))')

# Store a list of atoms near the G(N1) donors within a G(N1) dictionary.
g_n1_dict = {"don_index": [], "don_name": [], "don_resn": [], "don_resi": [], "don_chain": [],
             "acc_index": [], "acc_name": [], "acc_resn": [], "acc_resi": [], "acc_chain": []}
for donor in stored.g_n1_list:
    # Find the atoms near the donor.
    stored.nearby_atoms = []
    cmd.iterate(f'index {donor[0]} around {snakemake.config["search_dist"]}',
                'stored.nearby_atoms.append((index, name, resn, resi, chain))')
    # If the nearby atom is a C(N3), U(O2), or their deoxy counterparts, add them to the G(N1) dictionary along with the
    # info for the G(N1) donor.
    for atom in stored.nearby_atoms:
        if (atom[2], atom[1]) in [("C", "N3"), ("U", "O2"), ("DC", "N3"), ("DT", "O2")]:
            g_n1_dict["don_index"].append(donor[0])
            g_n1_dict["don_name"].append(donor[1])
            g_n1_dict["don_resn"].append(donor[2])
            g_n1_dict["don_resi"].append(donor[3])
            g_n1_dict["don_chain"].append(donor[4])
            g_n1_dict["acc_index"].append(atom[0])
            g_n1_dict["acc_name"].append(atom[1])
            g_n1_dict["acc_resn"].append(atom[2])
            g_n1_dict["acc_resi"].append(atom[3])
            g_n1_dict["acc_chain"].append(atom[4])

# Create a dataframe based on the G(N1) dictionary.
g_n1_df = pd.DataFrame(g_n1_dict)

# Acquire the H-bonding geometry measurements for the acceptors near each G(N1).
g_n1_h_bonds_dict = {"don_index": [], "don_name": [], "don_resn": [], "don_resi": [], "don_chain": [], "acc_index": [],
                     "acc_name": [], "acc_resn": [], "acc_resi": [], "acc_chain": [], "don_acc_distance": [],
                     "h_acc_distance": [], "h_angle": [], "h_dihedral": [], "h_name": []}
for atom_pair in g_n1_df.itertuples():
    # Store the donor and acceptor atom values for the dataframe row.
    don_list = [atom_pair.don_index, atom_pair.don_name, atom_pair.don_resn, atom_pair.don_resi, atom_pair.don_chain]
    acc_list = [atom_pair.acc_index, atom_pair.acc_name, atom_pair.acc_resn, atom_pair.acc_resi, atom_pair.acc_chain]
    # Retrieve the H-bond measurements for the atom pair.
    h_bond_list = eval_H_bonding.evaluate(don_list, acc_list, eq_class_mem_id, residue_library.RESIDUE_LIBRARY)
    # If the H-bond evaluation is not successful, print the error message(s) and exit.
    successful_completion = h_bond_list[0]
    if not successful_completion:
        for note in h_bond_list[1]:
            print(note)
        sys.exit(1)
    # Add the H-bond measurements to the dictionary.
    else:
        for h_bond in h_bond_list[1]:
            row = don_list + acc_list + h_bond
            idx = 0
            for key in g_n1_h_bonds_dict:
                g_n1_h_bonds_dict[key].append(row[idx])
                idx += 1
# Create a dataframe based on the dictionary.
g_n1_h_bonds_df = pd.DataFrame(g_n1_h_bonds_dict)

# Store a list of U(N3) donors (or its deoxy counterpart) from within or near the equivalence class member RNA chains.
stored.u_n3_list = []
cmd.iterate(f'({mem_rna_chains} expand {snakemake.config["search_dist"]}) and ((resn U or resn DT) and name N3)',
            'stored.u_n3_list.append((index, name, resn, resi, chain))')

# Store a list of atoms near the U(N3) donors within a U(N3) dictionary.
u_n3_dict = {"don_index": [], "don_name": [], "don_resn": [], "don_resi": [], "don_chain": [],
             "acc_index": [], "acc_name": [], "acc_resn": [], "acc_resi": [], "acc_chain": []}
for donor in stored.u_n3_list:
    # Find the atoms near the donor.
    stored.nearby_atoms = []
    cmd.iterate(f'index {donor[0]} around {snakemake.config["search_dist"]}',
                'stored.nearby_atoms.append((index, name, resn, resi, chain))')
    # If the nearby atom is a A(N1), G(O6), or their deoxy counterparts, add them to the U(N3) dictionary along with the
    # info for the U(N3) donor.
    for atom in stored.nearby_atoms:
        if (atom[2], atom[1]) in [("A", "N1"), ("G", "O6"), ("DA", "N1"), ("DG", "O6")]:
            u_n3_dict["don_index"].append(donor[0])
            u_n3_dict["don_name"].append(donor[1])
            u_n3_dict["don_resn"].append(donor[2])
            u_n3_dict["don_resi"].append(donor[3])
            u_n3_dict["don_chain"].append(donor[4])
            u_n3_dict["acc_index"].append(atom[0])
            u_n3_dict["acc_name"].append(atom[1])
            u_n3_dict["acc_resn"].append(atom[2])
            u_n3_dict["acc_resi"].append(atom[3])
            u_n3_dict["acc_chain"].append(atom[4])

# Create a dataframe based on the U(N3) dictionary.
u_n3_df = pd.DataFrame(u_n3_dict)

# Acquire the H-bonding geometry measurements for the acceptors near each U(N3).
u_n3_h_bonds_dict = {"don_index": [], "don_name": [], "don_resn": [], "don_resi": [], "don_chain": [], "acc_index": [],
                     "acc_name": [], "acc_resn": [], "acc_resi": [], "acc_chain": [], "don_acc_distance": [],
                     "h_acc_distance": [], "h_angle": [], "h_dihedral": [], "h_name": []}
for atom_pair in u_n3_df.itertuples():
    # Store the donor and acceptor atom values for the dataframe row.
    don_list = [atom_pair.don_index, atom_pair.don_name, atom_pair.don_resn, atom_pair.don_resi, atom_pair.don_chain]
    acc_list = [atom_pair.acc_index, atom_pair.acc_name, atom_pair.acc_resn, atom_pair.acc_resi, atom_pair.acc_chain]
    # Retrieve the H-bond measurements for the atom pair.
    h_bond_list = eval_H_bonding.evaluate(don_list, acc_list, eq_class_mem_id, residue_library.RESIDUE_LIBRARY)
    # If the H-bond evaluation is not successful, print the error message(s) and exit.
    successful_completion = h_bond_list[0]
    if not successful_completion:
        for note in h_bond_list[1]:
            print(note)
        sys.exit(1)
    # Add the H-bond measurements to the dictionary.
    else:
        for h_bond in h_bond_list[1]:
            row = don_list + acc_list + h_bond
            idx = 0
            for key in u_n3_h_bonds_dict:
                u_n3_h_bonds_dict[key].append(row[idx])
                idx += 1
# Create a dataframe based on the dictionary.
u_n3_h_bonds_df = pd.DataFrame(u_n3_h_bonds_dict)

# Create a dataframe of residues that are involved in potential H-bonding interactions.
res_df = pd.concat([
    nucleobase_df.loc[:, ["don_resn", "don_resi", "don_chain"]]
    .rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"}),
    nucleobase_df.loc[:, ["acc_resn", "acc_resi", "acc_chain"]]
    .rename(columns={"acc_resn": "resn", "acc_resi": "resi", "acc_chain": "chain"})]).drop_duplicates()

# Create a dataframe with the b-factors of heavy atoms for residues involved in potential H-bonding interactions.
b_factor_df = pd.DataFrame({"index": [], "name": [], "resn": [], "resi": [], "chain": [], "b-factor": [], "subset": []})
for res in res_df.itertuples():
    stored.atoms = []
    # Collect the b-factors for backbone atoms.
    cmd.iterate(f"resn {res.resn} and resi {res.resi} and chain {res.chain} and backbone and not elem H",
                "stored.atoms.append((index, name, resn, resi, chain, b, 'backbone'))")
    # Collect the b-factors for side chain atoms.
    cmd.iterate(f"resn {res.resn} and resi {res.resi} and chain {res.chain} and sidechain and not elem H",
                "stored.atoms.append((index, name, resn, resi, chain, b, 'sidechain'))")
    for atom in stored.atoms:
        b_factor_df.loc[len(b_factor_df)] = atom

# Record info on the atoms associated with the previously determined 100 atom indices and check whether there are any
# changes.
stored.check_two = []
for index in indices:
    cmd.iterate(f'index {index}', 'stored.check_two.append((index, name, resn, resi, chain))')
if not stored.check_one == stored.check_two:
    print(f"Error: The indices in equivalence class member {snakemake.wildcards.eq_class_members} changed.")
    sys.exit(1)

# Write a csv containing heavy atom to donor of interest distances.
with open(snakemake.output.distances, "w") as csv_file:
    writer = csv.writer(csv_file)
    if commit_hash:
        writer.writerow([f"# dual-H-bonding-nucleobases repo git commit hash: {commit_hash}"])
    writer.writerow([f"# input file: {snakemake.input[0]}"])
    writer.writerow([f"# file created on: {datetime.now().strftime('%y-%m-%d %H:%M:%S.%f')}"])
distance_df.to_csv(snakemake.output.distances, index=False, mode='a', na_rep='NaN')

# Write a csv containing all H-bonding information.
with open(snakemake.output.h_bond, "w") as csv_file:
    writer = csv.writer(csv_file)
    if commit_hash:
        writer.writerow([f"# dual-H-bonding-nucleobases repo git commit hash: {commit_hash}"])
    writer.writerow([f"# input file: {snakemake.input[0]}"])
    writer.writerow([f"# file created on: {datetime.now().strftime('%y-%m-%d %H:%M:%S.%f')}"])
pd.concat([don_h_bonds_df, g_n1_h_bonds_df, u_n3_h_bonds_df]).to_csv(snakemake.output.h_bond, index=False, mode='a',
                                                                     na_rep='NaN')

# Write a csv file that stores applicable nucleobases containing donors of interest.
nuc_grp = ["don_resn", "don_resi", "don_chain"]
with open(snakemake.output.nuc, "w") as csv_file:
    writer = csv.writer(csv_file)
    if commit_hash:
        writer.writerow([f"# dual-H-bonding-nucleobases repo git commit hash: {commit_hash}"])
    writer.writerow([f"# input file: {snakemake.input[0]}"])
    writer.writerow([f"# file created on: {datetime.now().strftime('%y-%m-%d %H:%M:%S.%f')}"])
nucleobase_df.loc[:, nuc_grp].drop_duplicates().to_csv(snakemake.output.nuc, index=False, mode='a', na_rep='NaN')

# Write a csv containing the b-factors.
with open(snakemake.output.b_factor, "w") as csv_file:
    writer = csv.writer(csv_file)
    if commit_hash:
        writer.writerow([f"# dual-H-bonding-nucleobases repo git commit hash: {commit_hash}"])
    writer.writerow([f"# input file: {snakemake.input[0]}"])
    writer.writerow([f"# file created on: {datetime.now().strftime('%y-%m-%d %H:%M:%S.%f')}"])
b_factor_df.to_csv(snakemake.output.b_factor, index=False, mode='a', na_rep='NaN')

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

# Save the modified structure if specified.
if snakemake.config["save_modified_mmcif"]:
    cmd.save(modified_mmcif_dir + eq_class_mem_id + ".cif")

# Remove the original mmCIF if specified.
if snakemake.config["remove_original_mmcif"]:
    subprocess.run(["rm", original_mmcif_dir + f"{eq_class_mem[1][0]}.cif".lower()])

# Close files and reset stdout and stderr.
stdout_file.close()
stderr_file.close()
sys.stdout = stdout
sys.stderr = stderr
