"""
This script iterates through the lines of a representative set file and collects the equivalence class names, PDB IDs,
model info, and chain info for the representative structures. The name of the representative set file must be provided
as the first argument, along with the folder where it resides (e.g., rep_set_files/nrlist_<identifier>.csv). The info
for each representative structure is written to a csv file in a folder named eq_class_files. If commit_hash is provided
as an argument, the commit hash of the repo will also be written within a commented line to the csv file if no
uncommitted changes have been made to the repo.
"""

import sys
import os
import csv
import subprocess
from datetime import datetime

# check whether an eq_class_files folder exists, and if it does not exist, create the folder
eq_class_dir = os.getcwd() + "/eq_class_files"
if not os.path.isdir(eq_class_dir):
    os.mkdir(eq_class_dir)

# If commit_hash is supplied as the first argument, check if any changes have been made to the repo and get the hash
# of the current git commit. If uncommitted changes have been made, print an error message and exit.
repo_changes = ""
commit_hash = ""
if "commit_hash" in sys.argv:
    repo_changes = subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=no"],
                                           cwd=os.path.dirname(os.path.realpath(__file__))).decode('ascii').strip()
    if not repo_changes:
        commit_hash = subprocess.check_output(["git", "rev-parse", "HEAD"],
                                              cwd=os.path.dirname(os.path.realpath(__file__))).decode('ascii').strip()
    else:
        print(f"Error: Uncommitted changes have been made to the repo.")
        sys.exit(1)

# check to make sure the specified representative set file exists
nrlist_file = ""
try:
    if not os.path.isfile(os.getcwd() + "/" + sys.argv[1]):
        print(f"Error: The file named {sys.argv[1]} was not found.")
        sys.exit(1)
    else:
        nrlist_file = os.getcwd() + "/" + sys.argv[1]
except IndexError:
    print(f"Error: The name of the representative set file must be provided as the first argument.")
    sys.exit(1)

# collect and write the information to csv files
with open(nrlist_file, mode='r') as read_file:
    for line in csv.reader(read_file):
        PDB_ID_list = []
        model_list = []
        chain_list = []
        for group in line[1].split('+'):
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
                print(f"Error: The representative structure for equivalence class {line[0]} is missing "
                      f"information.")
                sys.exit(1)
            # check if the PDB ID contains the expected number of characters
            if len("".join(PDB_ID)) != 4:
                print(f"Error: At least one of the PDB IDs for the representative structure of equivalence class "
                      f"{line[0]} does not contain the expected number of characters.")
                sys.exit(1)
            PDB_ID_list.append("".join(PDB_ID))
            model_list.append("".join(model))
            chain_list.append("".join(chain))
        # check to make sure the PDB IDs match
        for pdb in PDB_ID_list:
            if pdb != PDB_ID_list[0]:
                print(f"Error: The PDB IDs for the representative structure of equivalence class {line[0]} do not "
                      f"match.")
                sys.exit(1)
        with open(f"{eq_class_dir}/{line[0]}_info.csv", "w") as write_file:
            writer = csv.writer(write_file)
            if commit_hash:
                writer.writerow([f"# dual-H-bonding-nucleobases repo git commit hash: {commit_hash}"])
            writer.writerow([f"# input file: {nrlist_file}"])
            writer.writerow([f"# file created on: {datetime.now().strftime('%y-%m-%d %H:%M:%S.%f')}"])
            writer.writerow([line[0]])
            writer.writerow(PDB_ID_list)
            writer.writerow(model_list)
            writer.writerow(chain_list)
