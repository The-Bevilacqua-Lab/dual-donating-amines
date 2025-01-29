"""
Create a script that is to be used with PyMOL to review dual H-bonding amines from the combined.csv file that is output
from combine_data.py. This could also be a different CSV file with the same format as combined.csv. The dual H-bonding
amines in the data are then filtered based on the supplied arguments.
"""

import argparse
from datetime import datetime
import pandas as pd

# Parse the arguments.
parser = argparse.ArgumentParser(prog='Review Amines', description='Create a Python script that can be used with PyMOL '
                                                                   'to review identified dual H-bonding amines.')
parser.add_argument('input_file')
parser.add_argument('output_file')
parser.add_argument('--lower_eta', type=float,
                    help='What is the lower bound (inclusive) for the eta values, ranging from 0 to 360?')
parser.add_argument('--upper_eta', type=float,
                    help='What is the upper bound for the eta values, ranging from 0 to 360?')
parser.add_argument('--lower_theta', type=float,
                    help='What is the lower bound (inclusive) for the theta values, ranging from 0 to 360?')
parser.add_argument('--upper_theta', type=float,
                    help='What is the upper bound for the theta values, ranging from 0 to 360?')
parser.add_argument('--lower_chi', type=float,
                    help='What is the lower bound (inclusive) for the chi values?')
parser.add_argument('--upper_chi', type=float, help='What is the upper bound for the chi values?')
parser.add_argument('--don_resn', type=str,
                    help='What is the name of the donor residue (e.g., A, C, or G)?')
parser.add_argument('--search_dist', type=str,
                    help='What is the maximum distance from the parent nucleobase of the amine that surrounding '
                         'residues must reside within to be shown?')
args = parser.parse_args()

# Read the data.
df = pd.read_csv(args.input_file, comment="#", keep_default_na=False, na_values="NaN", dtype="str")

# Remove irrelevant or redundant rows such that each dual H-bonding amine is only represented one time.
df = df[df['type'] == '2'].drop_duplicates(subset=['don_index', 'eq_class_member'])

# Translate the pseudotorsions to range from 0 to 360 degrees.
df.loc[df['eta'].astype(float) >= 0, 'eta_translated'] = df.loc[df['eta'].astype(float) >= 0, 'eta'].astype(float)
df.loc[df['eta'].astype(float) < 0, 'eta_translated'] = df.loc[df['eta'].astype(float) < 0, 'eta'].astype(float) + 360
df.loc[df['theta'].astype(float) >= 0, 'theta_translated'] = (
    df.loc[df['theta'].astype(float) >= 0, 'theta'].astype(float))
df.loc[df['theta'].astype(float) < 0, 'theta_translated'] = (
        df.loc[df['theta'].astype(float) < 0, 'theta'].astype(float) + 360)

# Apply the filters specified by the arguments. Eta/theta values of -360 should not be included when eta/theta bounds
# are specified.
if args.lower_eta is not None:
    df = df[(df['eta'].astype(float) != -360) & (df['eta_translated'] >= args.lower_eta)]
if args.upper_eta is not None:
    df = df[(df['eta'].astype(float) != -360) & (df['eta_translated'] < args.upper_eta)]
if args.lower_theta is not None:
    df = df[(df['theta'].astype(float) != -360) & (df['theta_translated'] >= args.lower_theta)]
if args.upper_theta is not None:
    df = df[(df['theta'].astype(float) != -360) & (df['theta_translated'] < args.upper_theta)]
if args.lower_chi is not None:
    df = df[df['chi'].astype(float) >= args.lower_chi]
if args.upper_chi is not None:
    df = df[df['chi'].astype(float) < args.upper_chi]
if args.don_resn is not None:
    df = df[df['don_resn'] == args.don_resn]

# Write the Python script that is to be used with PyMOL.
with open(args.output_file, "w") as file:
    file.write(f"# file created on: {datetime.now().strftime('%y-%m-%d %H:%M:%S.%f')}")
    if args.lower_eta is not None:
        file.write(f"# eta lower bound (inclusive): {args.lower_eta}")
    if args.upper_eta is not None:
        file.write(f"# eta upper bound: {args.upper_eta}")
    if args.lower_theta is not None:
        file.write(f"# theta lower bound (inclusive): {args.lower_theta}")
    if args.upper_theta is not None:
        file.write(f"# theta upper bound: {args.upper_theta}")
    if args.lower_chi is not None:
        file.write(f"# chi lower bound (inclusive): {args.lower_chi}")
    if args.upper_chi is not None:
        file.write(f"# chi upper bound: {args.upper_chi}")
    if args.don_resn is not None:
        file.write(f"# don_resn: {args.don_resn}")
    file.write("# this script is to be used with PyMOL\n\n")
    file.write("from pymol import cmd\n\n")
    file.write("amines = [\n")
    for idx, row in enumerate(df.itertuples()):
        if idx == 0:
            file.write(f'    [[{row.don_name}, {row.don_resi}, {row.don_chain}, {row.acc_pair_1_name}, '
                       f'{row.acc_pair_1_resi}, {row.acc_pair_1_chain}, {row.PDB}, {row.model}]')
        else:
            file.write(f', [{row.don_name}, {row.don_resi}, {row.don_chain}, {row.acc_pair_2_name}, '
                       f'{row.acc_pair_2_resi}, {row.acc_pair_2_chain}, {row.PDB}, {row.model}]')
        if idx == df["don_index"].size - 1:
            file.write(']\n]\n\n')
        else:
            file.write('],\n')
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
    file.write("    cmd.show('sticks', 'byres all within " + {args.search_dist} +
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
