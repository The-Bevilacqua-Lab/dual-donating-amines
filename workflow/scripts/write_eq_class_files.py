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

# Redirect stdout and stderr to log files.
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

# Write the equivalence class name, PDB ID, model info, and chain info to a csv for the equivalence class member.
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
eq_class = "".join(eq_class_build)
pdb_id = "".join(pdb_id_build)
model = "".join(model_build)
chain_list = []
for chain in chain_build:
    chain_list.append("".join(chain))
with open(snakemake.output[0], "w") as write_file:
    writer = csv.writer(write_file)
    if commit_hash:
        writer.writerow([f"# dual-H-bonding-nucleobases repo git commit hash: {commit_hash}"])
    writer.writerow([f"# information based on: {snakemake.config['rep_set_file']}"])
    writer.writerow([f"# file created on: {datetime.now().strftime('%y-%m-%d %H:%M:%S.%f')}"])
    writer.writerow([eq_class])
    writer.writerow([pdb_id] * len(chain_list))
    writer.writerow([model] * len(chain_list))
    writer.writerow(chain_list)

# Close files and reset stdout and stderr.
stdout_file.close()
stderr_file.close()
sys.stdout = stdout
sys.stderr = stderr
