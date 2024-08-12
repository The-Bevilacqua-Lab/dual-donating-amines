"""
Build on the output from the collect_data.py script, adding columns to identify H-bonding interactions, assign each
donor of interest as no, single, or dual H-bonding, identify acceptors that bear a substantial negative charge, and
identify whether any given donor-acceptor pair will be used to inform the selection of H-bond geometry criteria.
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
    # If there are no H-bonding measurements, return NaN.
    if pd.isna(row["don_acc_distance"]) and pd.isna(row["h_acc_distance"]):
        return row
    # If the only distance measurement is between donor and acceptor, the donor is rotatable.
    elif pd.isna(row["h_acc_distance"]):
        if ((row["don_acc_distance"] <= DON_DIST_MAX) & (row["don_angle"] >= DON_ANG_MIN) &
                (row["don_angle"] <= DON_ANG_MAX)):
            row["h_bond"] = 1
            return row
        else:
            row["h_bond"] = 0
            return row
    # If both distance measurements exist, the donor is not rotatable.
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


# Read the data for the equivalence class member.
try:
    df = pd.read_csv(snakemake.input.data, comment="#", keep_default_na=False, na_values="NaN",
                     dtype={"don_resi": "object", "acc_resi": "object", "don_chain": "object", "acc_chain": "object"})
# Write an empty file and exit if one of the below errors is encountered.
except (FileNotFoundError, pd.errors.EmptyDataError):
    subprocess.run(["touch", snakemake.output.data])
    sys.exit(0)

# Set the H-bonding criteria.
H_DIST_MAX = snakemake.config["h_dist_max"]
H_ANG_MIN = snakemake.config["h_ang_min"]
DON_DIST_MAX = snakemake.config["don_dist_max"]
DON_ANG_MIN = snakemake.config["don_ang_min"]
DON_ANG_MAX = snakemake.config["don_ang_max"]

# Add a new column to identify whether each donor-acceptor pair meets the H-bond criteria.
df = df.assign(h_bond=pd.NA)
df = df.apply(h_bonding, axis=1)

# Add a new column that classifies the donors of interest as being no (0), single (1), or dual (2), H-bonding.
df = df.assign(type=pd.NA)
df.loc[:, ~df.columns.isin(["don_index"])] = (
    df.groupby("don_index", group_keys=False).apply(classify, include_groups=False))

# Prepare a list of acceptor residue and atom names from canonical amino acids and nucleic acids that have substantially
# greater negative charge.
can_neg = []
can_neut = []
for residue in residue_library.RESIDUE_LIBRARY:
    for acceptor in residue['acc']:
        if acceptor[3] < 0:
            can_neg.append(f"{residue['res']}.{acceptor[0]}")
        elif residue['res'] == 'ASP' and acceptor[0] == 'OD1':
            can_neg.append(f"{residue['res']}.{acceptor[0]}")
        elif residue['res'] == 'GLU' and acceptor[0] == 'OE1':
            can_neg.append(f"{residue['res']}.{acceptor[0]}")
        elif residue['res'] in ['A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT'] and acceptor[0] == 'OP1':
            can_neg.append(f"{residue['res']}.{acceptor[0]}")
        else:
            can_neut.append(f"{residue['res']}.{acceptor[0]}")

# Add a new column to identify acceptors from canonical amino acids and nucleic acids that have substantially greater
# negative charges.
df = df.assign(acc_charge=pd.NA)
df.loc[~df.isna()['acc_index'], "acc_charge"] = "other"
df.loc[df["acc_resn_name"].isin(can_neg), "acc_charge"] = "can_neg"
df.loc[df["acc_resn_name"].isin(can_neut), "acc_charge"] = "can_neut"
df = df.drop(columns=["don_resn_name", "acc_resn_name"])

# Add a new column to identify donor-acceptor pairs that will be used for evaluating H-bonding geometry. Each atom pair
# involves a donor of interest.
df = df.assign(geom=0)
don_acc_grp = ["don_index", "acc_index"]
(df.loc[
    # Only include atom pairs involving donors of interest.
    (df["DOI"] == 1) &
    # For a given donor-acceptor pair, only include the hydrogen with the smaller D-H...A distance.
    (df.groupby(don_acc_grp)["h_acc_distance"]
     .transform(lambda grp: [mem == grp.min() for mem in grp])) &
    # Do not consider A(N6)-U(O4), C(N4)-G(O6), G(N2)-C(O2), or G(N2)-C(N3) atom pairs.
    (~df[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["A", "N6", "U", "O4"])
     .all(axis='columns')) &
    (~df[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["C", "N4", "G", "O6"])
     .all(axis='columns')) &
    (~df[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["G", "N2", "C", "O2"])
     .all(axis='columns')) &
    (~df[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["G", "N2", "C", "N3"])
     .all(axis='columns')), "geom"]) = 1

# Write a csv containing the data with the new columns. Also report the H-bond criteria used.
with open(snakemake.output.data, "w") as csv_file:
    writer = csv.writer(csv_file)
    if commit_hash:
        writer.writerow([f"# dual-H-bonding-nucleobases repo git commit hash: {commit_hash}"])
    writer.writerow([f"# representative set file: {snakemake.config['rep_set_file']}"])
    writer.writerow([f"# file created on: {datetime.now().strftime('%y-%m-%d %H:%M:%S.%f')}"])
    writer.writerow([f"# H-acc max distance: {H_DIST_MAX}"])
    writer.writerow([f"# don-H-acc min angle: {H_ANG_MIN}"])
    writer.writerow([f"# don-acc max distance: {DON_DIST_MAX}"])
    writer.writerow([f"# don_ant-don-acc min angle: {DON_ANG_MIN}"])
    writer.writerow([f"# don_ant-don-acc max angle: {DON_ANG_MAX}"])
df.to_csv(snakemake.output.data, index=False, mode='a', na_rep='NaN')

# Close files and reset stdout and stderr.
stdout_file.close()
stderr_file.close()
sys.stdout = stdout
sys.stderr = stderr
