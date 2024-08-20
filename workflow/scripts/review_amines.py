"""
Create a script that is to be used with PyMOL to review dual H-bonding amines that involve a donor and acceptor that are
specified in the configuration file.
"""

import sys
import os
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


# Define a function to determine whether a donor of interest donates to an acceptor that is specified by the
# review_acceptor variable in the configuration file. If it does not, exclude the donor of interest from further
# consideration by returning an empty data frame.
def check_acceptor(group):
    nucleic_acids = ['A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT']
    keep_acc_resn = snakemake.config["review_acceptor"].split(".")[0]
    keep_acc_name = snakemake.config["review_acceptor"].split(".")[1]
    acc_present = False
    # Check whether any of the acceptors are the acceptor specified in the config file.
    for row in group.itertuples():
        # If the residue name is specified as "N", it can be any canonical nucleic acid.
        if ((row.acc_resn == keep_acc_resn or (keep_acc_resn == "N" and row.acc_resn in nucleic_acids)) and
                row.acc_name == keep_acc_name):
            acc_present = True
    # If the specified acceptor is not present, remove all rows in the group.
    if not acc_present:
        group = group.truncate(after=-1)
    return group


# Read the data from the different equivalence class members. Filter out structures and amines based on criteria
# specified within the configuration file.
for idx in range(len(snakemake.input.data)):
    try:
        df = pd.read_csv(snakemake.input.data[idx], comment="#", keep_default_na=False, na_values="NaN",
                         dtype={"don_resi": "object", "acc_resi": "object", "don_chain": "object",
                                "acc_chain": "object"})
        # If the data frame is empty, continue the loop.
        if df.empty:
            continue
        # Only keep H-bonding interactions for donors of interest that are dual H-bonding and that are specified by the
        # review_donor variable in the configuration file.
        keep_don_resn = snakemake.config["review_donor"].split(".")[0]
        keep_don_name = snakemake.config["review_donor"].split(".")[1]
        df = df[(df["h_bond"] == 1) & (df["DOI"] == 1) & (df["type"] == 2) & (df["don_resn"] == keep_don_resn) &
                (df["don_name"] == keep_don_name)]
        # Determine whether each donor of interest donates to an acceptor that is specified by the review_acceptor
        # variable in the configuration file.
        df.loc[:, ~df.columns.isin(["don_index"])] = (df.groupby("don_index", group_keys=False)
                                                      .apply(check_acceptor, include_groups=False))
        # Remove all rows that received NaN values from the check_acceptor function.
        df = df[~df["acc_index"].isna()]
        # Concatenate the data frame containing dual H-bonding amines from the single equivalence class member to the
        # data frame containing amines from all members.
        if idx == 0:
            combined_df = df.copy()
        else:
            # noinspection PyUnboundLocalVariable
            combined_df = pd.concat([combined_df, df])
    # Continue the loop if one of the below errors is encountered.
    except (FileNotFoundError, pd.errors.EmptyDataError):
        continue

# Write the Python script that is to be used with PyMOL.
with open(snakemake.output.script, "w") as file:
    if commit_hash:
        file.write(f"# dual-H-bonding-nucleobases repo git commit hash: {commit_hash}")
    file.write(f"# representative set file: {snakemake.config['rep_set_file']}")
    file.write(f"# file created on: {datetime.now().strftime('%y-%m-%d %H:%M:%S.%f')}")
    file.write("# this script is to be used with PyMOL\n\n")
    file.write("from pymol import cmd\n\n")
    file.write("amines = [\n")
    last_grp = list(combined_df.groupby(["don_index", "eq_class_member"]).groups.keys())[-1]
    for grp, mem in combined_df.groupby(["don_index", "eq_class_member"]):
        for i, pair in enumerate(mem.to_dict('index').values()):
            if i == 0:
                file.write(f'    [["{pair["don_name"]}", "{pair["don_resi"]}", "{pair["don_chain"]}", '
                           f'"{pair["acc_name"]}", "{pair["acc_resi"]}", "{pair["acc_chain"]}", "{pair["PDB"]}", '
                           f'"{pair["model"]}"]')
            else:
                file.write(f', ["{pair["don_name"]}", "{pair["don_resi"]}", "{pair["don_chain"]}", '
                           f'"{pair["acc_name"]}", "{pair["acc_resi"]}", "{pair["acc_chain"]}", "{pair["PDB"]}", '
                           f'"{pair["model"]}"]')
            if i == mem["don_index"].size - 1:
                if grp != last_grp:
                    file.write('],\n')
                else:
                    file.write(']\n]\n\n')
    file.write("# Establish initial settings.\n")
    file.write("amine_num = 1\n")
    file.write("current_pdb = '####'\n")
    file.write("current_model = -1\n\n\n")
    file.write("# Define a function to print the number of amines to review.\n")
    file.write("def get_num():\n")
    file.write("    print(f'Number of amines to review: {len(amines)}')\n\n\n")
    file.write("# Define a function to check if the current loaded PyMOL object matches the target object. If not, "
               "fetch the target\n")
    file.write("# structure. If the structure contains more than one model (or \"state\" in PyMOL), create a new "
               "object of just that\n")
    file.write("# model. If only a single model exists, change the name of the object to identify the model number.\n")
    file.write("def check_pdb_model(cur_pdb, cur_model, tar_pdb, tar_model):\n")
    file.write("    if not (cur_pdb == tar_pdb and cur_model == tar_model):\n")
    file.write("        cmd.delete('all')\n")
    file.write("        cmd.fetch(tar_pdb)\n")
    file.write("        if cmd.count_states(tar_pdb) > 1:\n")
    file.write("            cmd.create(f'{tar_pdb}_state_{tar_model}', selection=tar_pdb, source_state=tar_model,\n")
    file.write("                       target_state=1)\n")
    file.write("            cmd.delete(tar_pdb)\n")
    file.write("        else:\n")
    file.write("            cmd.set_name(tar_pdb, f'{tar_pdb}_state_{tar_model}')\n\n\n")
    file.write("# Define a function to advance to the amine matching the provided number designating its position in "
               "the amine list.\n")
    file.write("def goto(num):\n")
    file.write("    global current_pdb, current_model, amine_num\n")
    file.write("    num = int(num)\n")
    file.write("    amine_num = num\n")
    file.write("    if num < 1:\n")
    file.write("        print('The amine number must be 1 or greater.')\n")
    file.write("        return None\n")
    file.write("    i1 = num - 1\n")
    file.write("    target_pdb = amines[num-1][0][6]\n")
    file.write("    target_model = amines[num-1][0][7]\n")
    file.write("    check_pdb_model(current_pdb, current_model, target_pdb, target_model)\n")
    file.write("    current_pdb = target_pdb\n")
    file.write("    current_model = target_model\n")
    file.write("    cmd.delete(f'dist_{num}*')\n")
    file.write("    cmd.hide('all')\n")
    file.write("    cmd.color('grey60', 'elem C')\n")
    file.write("    cmd.show('sticks', f'resi {amines[num-1][0][1]} and chain {amines[num-1][0][2]}')\n")
    file.write("    cmd.color('cyan', 'elem C and visible')\n")
    file.write("    cmd.show('sticks', 'byres all within " + str(snakemake.config["search_dist"]) +
               " of (visible and sidechain)')\n")
    file.write("    for i2, pair in enumerate(amines[i1]):\n")
    file.write("        cmd.distance(f'dist_{num}_{i2}',\n"
               "                     f'name {amines[i1][i2][0]} and resi {amines[i1][i2][1]} "
               "and chain {amines[i1][i2][2]}',\n"
               "                     f'name {amines[i1][i2][3]} and resi {amines[i1][i2][4]} "
               "and chain {amines[i1][i2][5]}')\n")
    file.write("    cmd.hide('labels')\n")
    file.write("    cmd.show('dashes', f'dist_{num}*')\n")
    file.write("    cmd.hide('sticks', 'elem H')\n")
    file.write("    cmd.orient('visible')\n")
    file.write("    print(f'Showing amine number {num}')\n\n\n")
    file.write("# Show the next amine in the list.\n")
    file.write("def n():\n")
    file.write("    global amine_num\n")
    file.write("    if amine_num == len(amines):\n")
    file.write("        amine_num = 1\n")
    file.write("    else:\n")
    file.write("        amine_num += 1\n")
    file.write("    goto(amine_num)\n\n\n")
    file.write("# Show the previous amine in the list.\n")
    file.write("def p():\n")
    file.write("    global amine_num\n")
    file.write("    if amine_num == 1:\n")
    file.write("        amine_num = len(amines)\n")
    file.write("    else:\n")
    file.write("        amine_num -= 1\n")
    file.write("    goto(amine_num)\n\n\n")
    file.write("# Print the number of amines.\n")
    file.write("get_num()\n\n")
    file.write("# Load the first amine.\n")
    file.write("goto(amine_num)\n\n")
    file.write("# Extend the Python functions to act as PyMOL functions.\n")
    file.write("cmd.extend('get_num', get_num)\n")
    file.write("cmd.extend('goto', goto)\n")
    file.write("cmd.extend('n', n)\n")
    file.write("cmd.extend('p', p)\n")

# Close files and reset stdout and stderr.
stdout_file.close()
stderr_file.close()
sys.stdout = stdout
sys.stderr = stderr
