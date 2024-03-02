"""
This script iterates through the lines of a representative set file, searching for a specific equivalence class, and
collects the equivalence class name, PDB ID, model info, and chain info of the representative structures for the
specified equivalence class. The info for each representative structure is written to a csv file. The input file name,
output file name, and specified equivalence class are provided by the Snakefile. If commit_hash is set to true in the
Snakemake configuration file, the commit hash of the repo will also be written within a commented line to the csv file
if no uncommitted changes have been made to the repo.
"""

import sys
import os
import csv
import subprocess
from datetime import datetime

# redirect stdout and stderr to log files
stdout = sys.stdout
stderr = sys.stderr
stdout_file = open(snakemake.log.stdout, mode='w')
stderr_file = open(snakemake.log.stderr, mode='w')
sys.stdout = stdout_file
sys.stderr = stderr_file

# If commit_hash is supplied as the first argument, check if any changes have been made to the repo and get the hash
# of the current git commit. If uncommitted changes have been made to anything other than config/config.yaml, print an
# error message and exit.
repo_changes = []
commit_hash = ""
if snakemake.config["commit_hash"]:
    repo_changes = list(subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=no"],
                                                cwd=os.path.dirname(os.path.realpath(__file__)))
                        .decode('ascii').strip().split("\n"))
    if repo_changes == [""] or (len(repo_changes) == 1 and "config/config.yaml" in repo_changes[0]):
        commit_hash = subprocess.check_output(["git", "rev-parse", "HEAD"],
                                              cwd=os.path.dirname(os.path.realpath(__file__))).decode('ascii').strip()
    else:
        print(f"Error: Uncommitted changes have been made to the repo.")
        sys.exit(1)

# collect information for a specific equivalence class
eq_class = []
with open(snakemake.input[0], mode='r') as read_file:
    count = 0
    for line in csv.reader(read_file):
        if line[0] == snakemake.wildcards.eq_class:
            count += 1
            eq_class = line
    if count == 0:
        print(f"Error: The equivalence class {snakemake.wildcards.eq_class} was not found.")
        sys.exit(1)
    elif count > 1:
        print(f"Error: The equivalence class {snakemake.wildcards.eq_class} was found on more than one lines in the "
              f"representative set file.")
        sys.exit(1)

# parse the line describing the specified equivalence class
PDB_ID_list = []
model_list = []
chain_list = []
for group in eq_class[1].split('+'):
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
    # check to ensure this part of the representative structure contains the expected amount of information
    if index != 2:
        print(f"Error: The representative structure for equivalence class {eq_class[0]} is missing "
              f"information.")
        sys.exit(1)
    # check if the PDB ID contains the expected number of characters
    if len("".join(PDB_ID)) != 4:
        print(f"Error: At least one of the PDB IDs for the representative structure of equivalence class "
              f"{eq_class[0]} does not contain the expected number of characters.")
        sys.exit(1)
    PDB_ID_list.append("".join(PDB_ID))
    model_list.append("".join(model))
    chain_list.append("".join(chain))

# check to make sure the PDB IDs match
for pdb in PDB_ID_list:
    if pdb != PDB_ID_list[0]:
        print(f"Error: The PDB IDs for the representative structure of equivalence class {eq_class[0]} do not "
              f"match.")
        sys.exit(1)

# write the equivalence class name, PDB ID, model info, and chain info to a csv file
with open(snakemake.output[0], "w") as write_file:
    writer = csv.writer(write_file)
    if commit_hash:
        writer.writerow([f"# dual-H-bonding-nucleobases repo git commit hash: {commit_hash}"])
    writer.writerow([f"# input file: {snakemake.input[0]}"])
    writer.writerow([f"# file created on: {datetime.now().strftime('%y-%m-%d %H:%M:%S.%f')}"])
    writer.writerow([eq_class[0]])
    writer.writerow(PDB_ID_list)
    writer.writerow(model_list)
    writer.writerow(chain_list)

# close files and reset stdout and stderr
stdout_file.close()
stderr_file.close()
sys.stdout = stdout
sys.stderr = stderr
