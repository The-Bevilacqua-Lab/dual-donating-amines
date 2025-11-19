"""
Build on the output from the collect_data.py script, adding columns to identify H-bonding interactions, assign each
donor of interest as no, single, or dual H-bonding, and identify whether any given donor-acceptor pair will be used
to inform the selection of H-bond geometry criteria.
"""

import sys
import os
import subprocess
import csv
from datetime import datetime
import pandas as pd
import residue_library

# Redirect stdout and stderr to log files.
stdout = sys.stdout
stderr = sys.stderr
stdout_file = open(snakemake.log.stdout, mode='w')
stderr_file = open(snakemake.log.stderr, mode='w')
sys.stdout = stdout_file
sys.stderr = stderr_file


# This function is to be run if an error is encountered.
def error(msg):
    # Create the output files expected by Snakemake.
    subprocess.run(["touch", snakemake.output.data])
    # Print the error message.
    print(msg)
    # Close files, reset stdout and stderr, and exit.
    stdout_file.close()
    stderr_file.close()
    sys.stdout = stdout
    sys.stderr = stderr
    sys.exit(0)


# If commit_hash is set to true in the Snakemake configuration file, check if any changes have been made to the repo and
# get the hash of the current git commit. If uncommitted changes have been made to anything other than files listed in
# the acceptable_changes variable defined below, print an error message and exit.
repo_changes = []
commit_hash = ""
if snakemake.config["commit_hash"]:
    repo_changes = list(subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=no"],
                                                cwd=os.path.dirname(os.path.realpath(__file__)))
                        .decode('ascii').strip().split("\n"))
    acceptable_changes = ['config/config.yaml']
    for file in repo_changes:
        if file.split(' ')[-1] in acceptable_changes:
            repo_changes.remove(file)
    if len(repo_changes) == 0:
        commit_hash = subprocess.check_output(["git", "rev-parse", "HEAD"],
                                              cwd=os.path.dirname(os.path.realpath(__file__))).decode('ascii').strip()
    # Print the error message, close files, reset stdout and stderr, and exit.
    else:
        print("Error: Uncommitted changes have been made to the repo.")
        stdout_file.close()
        stderr_file.close()
        sys.stdout = stdout
        sys.stderr = stderr
        sys.exit(1)


# Define a function to determine whether a donor/acceptor atom pair meets the H-bonding criteria.
def h_bonding(row):
    # Return 0 for the h_bond column if no acceptor atom exists.
    if pd.isna(row["acc_resi"]):
        row["h_bond"] = 0
        return row
    # Issue an error message if the needed H-bonding measurements are NaN.
    elif pd.isna(row["h_acc_distance"]) or pd.isna(row["h_angle"]):
        error(f"Error: The needed H-bonding measurements are NaN for a donor/acceptor atom pair in {eq_class_mem_id}.")
    # If both measurements exist, evaluate them against the criteria.
    else:
        if (row["h_acc_distance"] <= H_DIST_MAX) & (row["h_angle"] >= H_ANG_MIN):
            row["h_bond"] = 1
            return row
        else:
            row["h_bond"] = 0
            return row


# Define a function to classify donors of interest as no (0), single (1), or dual (2), H-bonding. The donors that meet
# the dual classification must use two of its hydrogens in the H-bonding interactions which must involve at least two
# different acceptors.
def classify(grp):
    h_name_list = []
    acc_index_list = []
    for row in grp.itertuples():
        if not pd.isna(row.h_bond):
            if row.h_bond == 1 and row.h_name not in h_name_list:
                h_name_list.append(row.h_name)
            if row.h_bond == 1 and row.acc_index not in acc_index_list:
                acc_index_list.append(row.acc_index)
    if len(h_name_list) == 1 or len(acc_index_list) == 1:
        grp['type'] = 1
        return grp
    elif len(h_name_list) == 2 and len(acc_index_list) >= 2:
        grp['type'] = 2
        return grp
    else:
        grp['type'] = 0
        return grp


# Define a function to identify the pair of acceptors that a dual H-bonding amine is donating to. If there is more than
# one acceptor for a given hydrogen, choose the acceptor that exhibits the best (greatest) H-bonding angle.
def find_pairs(grp):
    # If the type is not 2, do not change the values.
    if pd.isna(grp.iloc[0, :].loc['type']) or grp.iloc[0, :].loc['type'] != 2:
        return grp
    h_name_list = []
    for row in grp.itertuples():
        if not pd.isna(row.h_bond):
            if row.h_bond == 1 and row.h_name not in h_name_list:
                h_name_list.append(row.h_name)
    if len(h_name_list) != 2:
        error("Error: During the process_data rule, exactly two hydrogens were not identified when the find_pairs "
              f"function was run for an amine classified as dual H-bonding in {eq_class_mem_id}.")
    acc_pair_1 = []
    acc_pair_2 = []
    for row in grp.itertuples():
        if row.h_bond == 1 and row.h_name == h_name_list[0]:
            acc_pair_1.append(row)
        if row.h_bond == 1 and row.h_name == h_name_list[1]:
            acc_pair_2.append(row)
    for pair in [acc_pair_1, acc_pair_2]:
        if len(pair) > 1:
            best_ang = pair[0]
            for named_series in pair[1:]:
                if named_series.h_angle > best_ang.h_angle:
                    best_ang = named_series
            if pair == acc_pair_1:
                acc_pair_1 = best_ang
            elif pair == acc_pair_2:
                acc_pair_2 = best_ang
        else:
            if pair == acc_pair_1:
                acc_pair_1 = pair[0]
            elif pair == acc_pair_2:
                acc_pair_2 = pair[0]
    # Do the acceptors belong to the same residue?
    same_resi = False
    if acc_pair_1.acc_resi == acc_pair_2.acc_resi and acc_pair_1.acc_chain == acc_pair_2.acc_chain:
        same_resi = True
    grp['acc_pair_1_name'], grp['acc_pair_1_resn'], grp['acc_pair_1_resi'], grp['acc_pair_1_chain'] = (
        acc_pair_1.acc_name, acc_pair_1.acc_resn, acc_pair_1.acc_resi, acc_pair_1.acc_chain)
    grp['acc_pair_2_name'], grp['acc_pair_2_resn'], grp['acc_pair_2_resi'], grp['acc_pair_2_chain'] = (
        acc_pair_2.acc_name, acc_pair_2.acc_resn, acc_pair_2.acc_resi, acc_pair_2.acc_chain)
    grp['same_resi'] = same_resi
    return grp


# Create a new column with the names of both atoms of the acceptor pair corresponding to dual H-bonding amines. Some
# names will be changed from what is used in the PDB. For instance, both OP1 and OP2 atoms will be renamed NPO. Which
# one of the two acceptors appears first will depend on how the two acceptors are sorted after combining into a list.
# Having consistent sorting of acceptors is important for the eventual plotting of this data. If the acceptor pairs
# are not consistently sorted, different instances of the same acceptor pair could be diluted into two bins. For
# instance, the pairings "A(N7), N(NPO)" and "N(NPO), A(N7)" would be binned separately when they are actually the same.
def combine_acc(grp):
    # If acc_pair_1_name is NA, do not change the values.
    if pd.isna(grp.iloc[0, :].loc['acc_pair_1_name']):
        return grp
    rna = ["A", "C", "G", "U"]
    dna = ["DA", "DC", "DG", "DT"]
    amino_acids = ["ALA", "ASP", "ASN", "ARG", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                   "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    # Create the name for the first acceptor.
    if (grp.iloc[0, :].loc['acc_pair_1_resn'] in rna and
            grp.iloc[0, :].loc['acc_pair_1_name'] in ["OP1", "OP2", "O2'", "O3'","O4'", "O5'"]):
        if grp.iloc[0, :].loc['acc_pair_1_name'] in ["OP1", "OP2"]:
            acc_pair_1_resn_name = "N(NPO)"
        else:
            acc_pair_1_resn_name = f"N({grp.iloc[0, :].loc['acc_pair_1_name']})"
    elif (grp.iloc[0, :].loc['acc_pair_1_resn'] in dna and
            grp.iloc[0, :].loc['acc_pair_1_name'] in ["OP1", "OP2", "O3'","O4'", "O5'"]):
        if grp.iloc[0, :].loc['acc_pair_1_name'] in ["OP1", "OP2"]:
            acc_pair_1_resn_name = "DN(NPO)"
        else:
            acc_pair_1_resn_name = f"DN({grp.iloc[0, :].loc['acc_pair_1_name']})"
    elif grp.iloc[0, :].loc['acc_pair_1_resn'] in amino_acids and grp.iloc[0, :].loc['acc_pair_1_name'] == "O":
        acc_pair_1_resn_name = "AA(O)"
    else:
        acc_pair_1_resn_name = f"{grp.iloc[0, :].loc['acc_pair_1_resn']}({grp.iloc[0, :].loc['acc_pair_1_name']})"
    # Create the name for the second acceptor.
    if (grp.iloc[0, :].loc['acc_pair_2_resn'] in rna and
            grp.iloc[0, :].loc['acc_pair_2_name'] in ["OP1", "OP2", "O2'", "O3'","O4'", "O5'"]):
        if grp.iloc[0, :].loc['acc_pair_2_name'] in ["OP1", "OP2"]:
            acc_pair_2_resn_name = "N(NPO)"
        else:
            acc_pair_2_resn_name = f"N({grp.iloc[0, :].loc['acc_pair_2_name']})"
    elif (grp.iloc[0, :].loc['acc_pair_2_resn'] in dna and
            grp.iloc[0, :].loc['acc_pair_2_name'] in ["OP1", "OP2", "O3'","O4'", "O5'"]):
        if grp.iloc[0, :].loc['acc_pair_2_name'] in ["OP1", "OP2"]:
            acc_pair_2_resn_name = "DN(NPO)"
        else:
            acc_pair_2_resn_name = f"DN({grp.iloc[0, :].loc['acc_pair_2_name']})"
    elif grp.iloc[0, :].loc['acc_pair_2_resn'] in amino_acids and grp.iloc[0, :].loc['acc_pair_2_name'] == "O":
        acc_pair_2_resn_name = "AA(O)"
    else:
        acc_pair_2_resn_name = f"{grp.iloc[0, :].loc['acc_pair_2_resn']}({grp.iloc[0, :].loc['acc_pair_2_name']})"
    # Combine the names within a list, then sort it.
    combined_list = [acc_pair_1_resn_name, acc_pair_2_resn_name]
    combined_list.sort()
    # Update the acc_pair_combined column to a string representation of the sorted list and then return.
    grp['acc_pair_combined'] = f"{combined_list[0]}, {combined_list[1]}"
    return grp

# Create an identifier for the equivalence class member.
eq_class_mem_id = snakemake.input.data[21:-4]

# Read the data for the equivalence class member.
try:
    df = pd.read_csv(snakemake.input.data, comment="#", keep_default_na=False, na_values="NaN",
                     dtype={"don_resi": "str", "acc_resi": "str", "don_chain": "str", "acc_chain": "str", "PDB": "str",
                            "PDB_version_number": "str", "PDB_version_date": "str"})
# Write an empty file and exit if FileNotFoundError is encountered.
except FileNotFoundError:
    subprocess.run(["touch", snakemake.output.data])
    print(f"Warning: A FileNotFoundError was encountered with {eq_class_mem_id}.")
    sys.exit(0)
# Write an empty file and exit if EmptyDataError is encountered.
except pd.errors.EmptyDataError:
    subprocess.run(["touch", snakemake.output.data])
    print(f"Note: The collect_data output file for {eq_class_mem_id} contains no data.")
    sys.exit(0)

# Set the H-bonding criteria.
H_DIST_MAX = snakemake.config["h_dist_max"]
H_ANG_MIN = snakemake.config["h_ang_min"]

# Add a new column to identify whether each donor-acceptor pair meets the H-bond criteria.
df = df.assign(h_bond=pd.NA)
df = df.apply(h_bonding, axis=1)

# Add a new column that classifies the donors of interest as being no (0), single (1), or dual (2), H-bonding.
df = df.assign(type=pd.NA)
df.loc[:, ~df.columns.isin(["don_index"])] = (
    df.groupby("don_index", group_keys=False).apply(classify, include_groups=False))

# Add new columns that identifies the pair of acceptors corresponding to dual H-bonding amines.
df = df.assign(acc_pair_1_name=pd.NA, acc_pair_1_resn=pd.NA, acc_pair_1_resi=pd.NA, acc_pair_1_chain=pd.NA,
               acc_pair_2_name=pd.NA, acc_pair_2_resn=pd.NA, acc_pair_2_resi=pd.NA, acc_pair_2_chain=pd.NA,
               same_resi=pd.NA)
df.loc[:, ~df.columns.isin(["don_index"])] = (
    df.groupby("don_index", group_keys=False).apply(find_pairs, include_groups=False))

# Create a new column with the names of both acceptors. See the combine_acc function for more details.
df = df.assign(acc_pair_combined=pd.NA)
df.loc[:, ~df.columns.isin(["don_index"])] = (
    df.groupby("don_index", group_keys=False).apply(combine_acc, include_groups=False))

# Add a new column to identify donor-acceptor pairs that will be used for evaluating H-bonding geometry. Each atom pair
# involves a donor of interest.
df = df.assign(geom=0)
don_acc_grp = ["don_index", "acc_index"]
(df.loc[
    # For a given donor-acceptor pair, only include the hydrogen with the greater D-H...A angle.
    (df.groupby(don_acc_grp)["h_angle"]
     .transform(lambda grp: [mem == grp.max() for mem in grp])) &
    # Do not consider A(N6)-U(O4), C(N4)-G(O6), G(N2)-C(O2), or G(N2)-C(N3) atom pairs.
    (~df[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["A", "N6", "U", "O4"])
     .all(axis='columns')) &
    (~df[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["C", "N4", "G", "O6"])
     .all(axis='columns')) &
    (~df[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["G", "N2", "C", "O2"])
     .all(axis='columns')) &
    (~df[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["G", "N2", "C", "N3"])
     .all(axis='columns')), "geom"]) = 1

# Use the previously added column to identify the donor-acceptor pairs that were excluded.
(df.loc[
    # For a given donor-acceptor pair, only include the hydrogen with the greater D-H...A angle.
    (df.groupby(don_acc_grp)["h_angle"]
     .transform(lambda grp: [mem == grp.max() for mem in grp])) &
    # Only include A(N6)-U(O4) atom pairs.
    (df[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["A", "N6", "U", "O4"])
     .all(axis='columns')), "geom"]) = 2
(df.loc[
    # For a given donor-acceptor pair, only include the hydrogen with the greater D-H...A angle.
    (df.groupby(don_acc_grp)["h_angle"]
     .transform(lambda grp: [mem == grp.max() for mem in grp])) &
    # Only include C(N4)-G(O6) atom pairs.
    (df[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["C", "N4", "G", "O6"])
     .all(axis='columns')), "geom"]) = 3
(df.loc[
    # For a given donor-acceptor pair, only include the hydrogen with the greater D-H...A angle.
    (df.groupby(don_acc_grp)["h_angle"]
     .transform(lambda grp: [mem == grp.max() for mem in grp])) &
    # Only include G(N2)-C(O2) atom pairs.
    (df[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["G", "N2", "C", "O2"])
     .all(axis='columns')), "geom"]) = 4
(df.loc[
    # For a given donor-acceptor pair, only include the hydrogen with the greater D-H...A angle.
    (df.groupby(don_acc_grp)["h_angle"]
     .transform(lambda grp: [mem == grp.max() for mem in grp])) &
    # Only include G(N2)-C(N3) atom pairs.
    (df[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["G", "N2", "C", "N3"])
     .all(axis='columns')), "geom"]) = 5

# Remove columns that do not need to be written to the csv.
df = df.drop(columns=["don_resn_name", "acc_resn_name"])

# Write a csv containing the data with the new columns. Also report the H-bond criteria used.
with open(snakemake.output.data, "w") as csv_file:
    writer = csv.writer(csv_file)
    if commit_hash:
        writer.writerow([f"# dual-donating-amines repo git commit hash: {commit_hash}"])
    writer.writerow([f"# representative set file: {snakemake.config['ifes_file']}"])
    writer.writerow([f"# file created on: {datetime.now().strftime('%y-%m-%d %H:%M:%S.%f')}"])
    writer.writerow([f"# H-acc max distance: {H_DIST_MAX}"])
    writer.writerow([f"# don-H-acc min angle: {H_ANG_MIN}"])
df.to_csv(snakemake.output.data, index=False, mode='a', na_rep='NaN')

# Close files and reset stdout and stderr.
stdout_file.close()
stderr_file.close()
sys.stdout = stdout
sys.stderr = stderr
