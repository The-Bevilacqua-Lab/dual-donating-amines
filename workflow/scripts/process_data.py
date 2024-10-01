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

# TODO incorporate b-factor filter including removing residues with b-factor of 0

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

# Create a data frame of H-bonding atom pairs that include an acceptor of interest.
aoi_h_bond_df = df[(df["acc_resn_name"].isin(snakemake.config["acceptors_of_interest"])) & (df["h_bond"] == 1)].copy()

# Create a data frame of H-bonding atom pairs involving a dual H-bonding amine and where the nucleobase of the amine
# does not form H-bonds via any of its acceptors of interest. These are candidates for QM calculations.
qm_df = (df[(df["don_resn_name"].isin(snakemake.config["donors_of_interest"])) &
            (df["h_bond"] == 1) & (df["type"] == 2)]
         .merge(aoi_h_bond_df, indicator=True, how='left',
                left_on=['don_resn', 'don_resi', 'don_chain'],
                right_on=['acc_resn', 'acc_resi', 'acc_chain']))
qm_df = (qm_df[qm_df["_merge"] == "left_only"].assign(qm=1).rename(columns={"don_index_x": "don_index"})
         .loc[:, ["don_index", "qm"]]).drop_duplicates()
if not qm_df.empty:
    df = df.merge(qm_df, how='outer')
else:
    df = df.assign(qm=pd.NA)

# Create data frames of specific H-bonding interactions.
n6_o4_h_bond_df = (df[
                       (df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                        .eq(["N6", "A", "O4", "U", 1]).all(axis='columns')) |
                       (df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                        .eq(["N6", "DA", "O4", "U", 1]).all(axis='columns')) |
                       (df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                        .eq(["N6", "A", "O4", "DT", 1]).all(axis='columns')) |
                       (df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                        .eq(["N6", "DA", "O4", "DT", 1]).all(axis='columns'))
                   ].rename(columns=lambda col: f'{col}_N6_O4')).copy()
n3_n1_h_bond_df = (df[
                       (df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                        .eq(["N3", "U", "N1", "A", 1]).all(axis='columns')) |
                       (df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                        .eq(["N3", "DT", "N1", "A", 1]).all(axis='columns')) |
                       (df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                        .eq(["N3", "U", "N1", "DA", 1]).all(axis='columns')) |
                       (df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                        .eq(["N3", "DT", "N1", "DA", 1]).all(axis='columns'))
                   ].rename(columns=lambda col: f'{col}_N3_N1')).copy()
n4_o6_h_bond_df = (df[
                       (df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                        .eq(["N4", "C", "O6", "G", 1]).all(axis='columns')) |
                       (df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                        .eq(["N4", "DC", "O6", "G", 1]).all(axis='columns')) |
                       (df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                        .eq(["N4", "C", "O6", "DG", 1]).all(axis='columns')) |
                       (df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                        .eq(["N4", "DC", "O6", "DG", 1]).all(axis='columns'))
                   ].rename(columns=lambda col: f'{col}_N4_O6')).copy()
n1_n3_h_bond_df = (df[
                       (df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                        .eq(["N1", "G", "N3", "C", 1]).all(axis='columns')) |
                       (df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                        .eq(["N1", "DG", "N3", "C", 1]).all(axis='columns')) |
                       (df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                        .eq(["N1", "G", "N3", "DC", 1]).all(axis='columns')) |
                       (df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                        .eq(["N1", "DG", "N3", "DC", 1]).all(axis='columns'))
                   ].rename(columns=lambda col: f'{col}_N1_N3')).copy()
n2_o2_h_bond_df = (df[
                       (df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                        .eq(["N2", "G", "O2", "C", 1]).all(axis='columns')) |
                       (df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                        .eq(["N2", "DG", "O2", "C", 1]).all(axis='columns')) |
                       (df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                        .eq(["N2", "G", "O2", "DC", 1]).all(axis='columns')) |
                       (df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                        .eq(["N2", "DG", "O2", "DC", 1]).all(axis='columns'))
                   ].rename(columns=lambda col: f'{col}_N2_O2')).copy()

# Prepare a dataframe of nucleobases containing a donor of interest and involved in a canonical base pair.
base_pair_df = pd.DataFrame(columns=["don_index", "base_pair"])
# AU base pair
if not (n6_o4_h_bond_df.empty or n3_n1_h_bond_df.empty):
    base_pair_df = pd.concat([base_pair_df,
                              (n6_o4_h_bond_df[n6_o4_h_bond_df["DOI_N6_O4"] == 1]
                               .merge(n3_n1_h_bond_df,
                                      left_on=["don_resn_N6_O4", "don_resi_N6_O4", "don_chain_N6_O4", "acc_resn_N6_O4",
                                               "acc_resi_N6_O4", "acc_chain_N6_O4"],
                                      right_on=["acc_resn_N3_N1", "acc_resi_N3_N1", "acc_chain_N3_N1", "don_resn_N3_N1",
                                                "don_resi_N3_N1", "don_chain_N3_N1"], how='inner')
                               .rename(columns={"don_index_N6_O4": "don_index"})
                               .assign(base_pair="AU").loc[:, ["don_index", "base_pair"]])])
# CG base pair
if not (n4_o6_h_bond_df.empty or n1_n3_h_bond_df.empty or n2_o2_h_bond_df.empty):
    base_pair_df = pd.concat([base_pair_df,
                              (n4_o6_h_bond_df[n4_o6_h_bond_df["DOI_N4_O6"] == 1]
                               .merge(n1_n3_h_bond_df,
                                      left_on=["don_resn_N4_O6", "don_resi_N4_O6", "don_chain_N4_O6", "acc_resn_N4_O6",
                                               "acc_resi_N4_O6", "acc_chain_N4_O6"],
                                      right_on=["acc_resn_N1_N3", "acc_resi_N1_N3", "acc_chain_N1_N3", "don_resn_N1_N3",
                                                "don_resi_N1_N3", "don_chain_N1_N3"], how='inner')
                               .merge(n2_o2_h_bond_df,
                                      left_on=["don_resn_N4_O6", "don_resi_N4_O6", "don_chain_N4_O6", "acc_resn_N4_O6",
                                               "acc_resi_N4_O6", "acc_chain_N4_O6"],
                                      right_on=["acc_resn_N2_O2", "acc_resi_N2_O2", "acc_chain_N2_O2", "don_resn_N2_O2",
                                                "don_resi_N2_O2", "don_chain_N2_O2"], how='inner')
                               .rename(columns={"don_index_N4_O6": "don_index"})
                               .assign(base_pair="CG").loc[:, ["don_index", "base_pair"]])])
# GC base pair
if not (n4_o6_h_bond_df.empty or n1_n3_h_bond_df.empty or n2_o2_h_bond_df.empty):
    base_pair_df = pd.concat([base_pair_df,
                              (n2_o2_h_bond_df[n2_o2_h_bond_df["DOI_N2_O2"] == 1]
                               .merge(n1_n3_h_bond_df,
                                      left_on=["don_resn_N2_O2", "don_resi_N2_O2", "don_chain_N2_O2", "acc_resn_N2_O2",
                                               "acc_resi_N2_O2", "acc_chain_N2_O2"],
                                      right_on=["don_resn_N1_N3", "don_resi_N1_N3", "don_chain_N1_N3", "acc_resn_N1_N3",
                                                "acc_resi_N1_N3", "acc_chain_N1_N3"], how='inner')
                               .merge(n4_o6_h_bond_df,
                                      left_on=["don_resn_N2_O2", "don_resi_N2_O2", "don_chain_N2_O2", "acc_resn_N2_O2",
                                               "acc_resi_N2_O2", "acc_chain_N2_O2"],
                                      right_on=["acc_resn_N4_O6", "acc_resi_N4_O6", "acc_chain_N4_O6", "don_resn_N4_O6",
                                                "don_resi_N4_O6", "don_chain_N4_O6"], how='inner')
                               .rename(columns={"don_index_N2_O2": "don_index"})
                               .assign(base_pair="GC").loc[:, ["don_index", "base_pair"]])])
df = df.merge(base_pair_df, how='outer')

# Remove columns that do not need to be written to the csv.
df = df.drop(columns=["don_resn_name", "acc_resn_name"])

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
