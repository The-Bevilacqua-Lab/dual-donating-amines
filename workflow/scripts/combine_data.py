"""
Combine all the output data files from the process_data.py script.
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
    acceptable_changes = ['config/config.yaml', snakemake.config["rep_set_file"]]
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
combined_df = pd.DataFrame(columns=['don_index', 'don_name', 'don_resn', 'don_resi', 'don_chain', 'don_segi', 'don_alt',
                                    'count_1', 'count_2', 'b_factor', 'chi', 'acc_index', 'acc_name', 'acc_resn',
                                    'acc_resi', 'acc_chain', 'acc_segi', 'acc_alt', 'don_acc_distance',
                                    'h_acc_distance', 'don_angle', 'acc_angle', 'h_angle', 'h_dihedral', 'h_name',
                                    'model', 'PDB', 'eq_class_member', 'h_bond', 'type', 'acc_charge', 'geom', 'v0',
                                    'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'eta', 'theta', "eta'", "theta'", 'alpha',
                                    'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'i1'])

# Combine the dataframes from the different equivalence class members. Do not incorporate empty dataframes.
for idx in range(len(snakemake.input.data)):
    try:
        member_df = pd.read_csv(snakemake.input.data[idx], comment="#", keep_default_na=False, na_values="NaN",
                                dtype={"don_resi": "str", "acc_resi": "str", "don_chain": "str",
                                       "acc_chain": "str"})
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

# Close files and reset stdout and stderr.
stdout_file.close()
stderr_file.close()
sys.stdout = stdout
sys.stderr = stderr
