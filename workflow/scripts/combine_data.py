"""
Combine all the output data files from the process_data.py script and perform some additional operations to enable
analysis of donors of interest involved in canonical base pairing and H-bond distances to acceptors of interest.
"""

import sys
import os
import csv
import subprocess
from datetime import datetime
import pandas as pd

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

# Initialize the combined dataframe.
combined_df = pd.DataFrame(columns=['don_index', 'don_name', 'don_resn', 'don_resi', 'don_chain', 'don_segi', 'count_1',
                                    'count_2', 'b_factor', 'DOI', 'don_can_NA', 'acc_index', 'acc_name', 'acc_resn',
                                    'acc_resi', 'acc_chain', 'acc_segi', 'don_acc_distance', 'h_acc_distance',
                                    'don_angle', 'h_angle', 'h_dihedral', 'h_name', 'AOI', 'model', 'PDB',
                                    'eq_class_member', 'h_bond', 'type', 'acc_charge', 'geom', 'base_pair'])

# Combine the dataframes from the different equivalence class members. Do not incorporate empty dataframes.
for idx in range(len(snakemake.input.data)):
    try:
        member_df = pd.read_csv(snakemake.input.data[idx], comment="#", keep_default_na=False, na_values="NaN",
                                dtype={"don_resi": "object", "acc_resi": "object", "don_chain": "object",
                                       "acc_chain": "object"})
        if combined_df.empty and not member_df.empty:
            combined_df = member_df.copy()
        elif not member_df.empty:
            combined_df = pd.concat([combined_df, member_df])
    except (FileNotFoundError, pd.errors.EmptyDataError):
        continue
combined_df = combined_df.reset_index(drop=True)

# Write the combined data to a csv file.
with open(snakemake.output.combined, "w") as csv_file:
    writer = csv.writer(csv_file)
    if commit_hash:
        writer.writerow([f"# dual-H-bonding-nucleobases repo git commit hash: {commit_hash}"])
    writer.writerow([f"# representative set file: {snakemake.config['rep_set_file']}"])
    writer.writerow([f"# file created on: {datetime.now().strftime('%y-%m-%d %H:%M:%S.%f')}"])
combined_df.to_csv(snakemake.output.combined, index=False, mode='a', na_rep='NaN')


# Define a function to determine whether the residue donating to the acceptor of interest is the same as the residue
# accepting from the donor of interest. If it is the same residue, the interactions are considered to be overlapping.
# The acceptor of interest and donor of interest are a part of the same nucleobase. Additionally, identify whether any
# of the atoms accepting an H-bond from the donor of interest bear substantial negative charge.
def check_overlap_charge(grp):
    overlap = 0
    neg_count = 0
    if any(grp["type"] == 0):
        grp.loc[:, ['overlap', "neg_charge"]] = [0, pd.NA]
        return grp
    for row in grp.itertuples():
        # Is the residue accepting the H-bond from the amine the same residue that is donating an H-bond to the acceptor
        # of interest?
        if row.acc_resn == row.don_resn_AOI and row.acc_resi == row.don_resi_AOI and row.acc_chain == row.don_chain_AOI:
            # Are the acceptor and donor atoms both part of the nucleobase (e.i., not backbone atoms)?
            if row.acc_name not in ["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"] and row.don_name_AOI not in ["O2'"]:
                overlap = 1
        if row.acc_charge == "can_neg":
            neg_count += 1
    grp['overlap'] = overlap
    # The DOI does not donate to any anionic acceptors.
    if neg_count == 0:
        grp["neg_charge"] = "None"
    # The DOI donates to neutral and anionic acceptors.
    elif neg_count < len(grp):
        grp["neg_charge"] = "Partial"
    # The DOI only donates to anionic acceptors.
    elif neg_count == len(grp):
        grp["neg_charge"] = "All"
    return grp


# Define a function to determine whether and donor of interest donates to any residues not included in the
# included_residues list specified within the config file. If it does, exclude the donor of interest from further
# consideration by returning an empty data frame.
def check_included_res(grp):
    for row in grp.itertuples():
        if row.acc_resn not in snakemake.config["included_residues"]:
            grp = grp.truncate(after=-1)
            return grp
    return grp


# Create a data frame that includes H-bonding atom pairs that include an acceptor of interest.
aoi_df = (combined_df[(combined_df["AOI"] == 1) & (combined_df["h_bond"] == 1)]
          .drop_duplicates(subset=["acc_index", "don_index", "eq_class_member"])
          .loc[:, ['don_index', 'don_name', 'don_resn', 'don_resi', 'don_chain', 'acc_index', 'acc_name', 'acc_resn',
                   'acc_resi', 'acc_chain', 'eq_class_member', 'don_acc_distance']]
          .copy())

# Create a data frame that includes H-bonding atom pairs that include a donor of interest which does not donate to any
# residues, or only those included in the included_residues list specified in the config file.
doi_df = combined_df[(combined_df["DOI"] == 1) & (combined_df["type"] > 0)].copy()
doi_df.loc[:, ~doi_df.columns.isin(["don_index", "eq_class_member"])] = (
    doi_df.groupby(['don_index', 'eq_class_member'], group_keys=False).apply(check_included_res, include_groups=False))
doi_df = pd.concat([doi_df, combined_df[(combined_df["DOI"] == 1) & (combined_df["type"] == 0)]])

# Create a merged data frame where residues bearing donors of interest are matched with residues bearing acceptors of
# interest.
merged_df = (doi_df.merge(aoi_df, left_on=['don_resn', 'don_resi', 'don_chain', 'eq_class_member'],
                          right_on=['acc_resn', 'acc_resi', 'acc_chain', 'eq_class_member'], suffixes=[None, "_AOI"]))

# Add overlap and neg_charge columns to the merged data frame. Refer to the check_overlap_charge function definition
# above for what these columns indicate.
merged_df.loc[:, ["overlap", "neg_charge"]] = pd.NA
merged_df.loc[:, ~merged_df.columns.isin(['don_index_AOI', "acc_index_AOI", "eq_class_member"])] = (
    merged_df.groupby(['don_index_AOI', 'acc_index_AOI', 'eq_class_member'], group_keys=False)
    .apply(check_overlap_charge, include_groups=False))

# Remove unnecessary columns and redundant rows from the merged data frame.
drop_columns = ["don_index", "don_resi", "don_chain", "don_segi", "count_1", "count_2", "DOI", "don_can_NA",
                "acc_index", "acc_name", "acc_resn", "acc_resi", "acc_chain", "acc_segi", "don_acc_distance",
                "h_acc_distance", "don_angle", "h_angle", "h_dihedral", "h_name", "AOI", "h_dihedral", "model", "PDB",
                "eq_class_member", "acc_charge", "geom", "h_bond", "don_index_AOI", "don_resi_AOI", "don_chain_AOI",
                "acc_index_AOI", "acc_resn_AOI", "acc_resi_AOI", "acc_chain_AOI"]
merged_df = (merged_df.drop_duplicates(subset=['don_index_AOI', 'acc_index_AOI', 'eq_class_member'])
             .drop(columns=drop_columns))

# Write the merged data frame to a csv file.
with open(snakemake.output.distances, "w") as csv_file:
    writer = csv.writer(csv_file)
    if commit_hash:
        writer.writerow([f"# dual-H-bonding-nucleobases repo git commit hash: {commit_hash}"])
    writer.writerow([f"# representative set file: {snakemake.config['rep_set_file']}"])
    writer.writerow([f"# file created on: {datetime.now().strftime('%y-%m-%d %H:%M:%S.%f')}"])
merged_df.to_csv(snakemake.output.distances, index=False, mode='a', na_rep='NaN')

# Close files and reset stdout and stderr.
stdout_file.close()
stderr_file.close()
sys.stdout = stdout
sys.stderr = stderr
