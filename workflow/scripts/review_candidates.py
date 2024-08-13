"""
Create a script that is to be used with PyMOL to review identified candidates.
"""

import subprocess
import sys
import pandas as pd

# Redirect stdout and stderr to log files.
stdout = sys.stdout
stderr = sys.stderr
stdout_file = open(snakemake.log.stdout, mode='w')
stderr_file = open(snakemake.log.stderr, mode='w')
sys.stdout = stdout_file
sys.stderr = stderr_file

# Read the data on the candidates.
df = pd.read_csv(snakemake.input.candidates, comment="#", keep_default_na=False, na_values="NaN",
                 dtype={"don_resi": "object", "don_chain": "object"})

# Create a list of candidates to cycle through.
candidates = []
for row in df.itertuples():
    candidates.append([row.PDB, row.model, row.don_resi, row.don_chain, row.b_factor])

# Write the Python script that is to be used with PyMOL.
with open(snakemake.output.script, "w") as file:
    file.write("# this script is to be used with PyMOL\n")
    file.write(f"# info on interactions came from {snakemake.input.candidates}\n")
    file.write("from pymol import cmd\n")
    file.write("from pymol import stored\n\n")
    file.write("candidates = [\n")
    for i, candidate in enumerate(candidates):
        file.write(f'    {candidate}')
        if i != len(candidates) - 1:
            file.write(',\n')
        else:
            file.write('\n]\n\n')
    file.write("# Establish initial settings.\n")
    file.write("cmd.set('scene_animation_duration', '0.5')\n")
    file.write("cmd.set('label_size', '-0.7')\n")
    file.write("cmd.set('label_font_id', '7')\n")
    file.write("cmd.set('label_color', 'white')\n")
    file.write("cmd.set('label_outline_color', 'black')\n")
    file.write("can_num = 1\n")
    file.write("current_pdb = '####'\n")
    file.write("current_model = -1\n\n\n")
    file.write("# Print the number of candidates.\n")
    file.write("def get_num():\n")
    file.write("    print(f'Number of candidates: {len(candidates)}')\n\n\n")
    file.write("# Check if the current loaded PyMOL object matches the target object. If not, fetch the target "
               "structure. If the\n")
    file.write("# structure contains more than one model (or \"state\" in PyMOL), create a new object of just that "
               "model. If only a single\n")
    file.write("# model exists, change the name of the object to identify the model number.\n")
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
    file.write("# Label the b-factor at the middle of the 6-membered ring.\n")
    file.write("def label_b_factor(res, b_factor):\n")
    file.write("    stored.suffix = ' \u00c5\u00b2'\n")
    file.write("    cmd.delete('b_factor_loc')\n")
    file.write("    cmd.pseudoatom('b_factor_loc', f'({res}) and (name N1 name C2 name N3 name C4 name C5 name C6)')\n")
    file.write("    cmd.hide('everything', 'b_factor_loc')\n")
    file.write("    cmd.label('b_factor_loc', f'str(round({b_factor}, 1)) + stored.suffix')\n\n\n")
    file.write("# Advance to the candidate matching the provided number designating its position in the candidate "
               "list.\n")
    file.write("def goto(num):\n")
    file.write("    global current_pdb, current_model\n")
    file.write("    num = int(num)\n")
    file.write("    if num < 1:\n")
    file.write("        print('The candidate number must be 1 or greater.')\n")
    file.write("        return None\n")
    file.write("    target_pdb = candidates[num-1][0]\n")
    file.write("    target_model = candidates[num-1][1]\n")
    file.write("    check_pdb_model(current_pdb, current_model, target_pdb, target_model)\n")
    file.write("    current_pdb = target_pdb\n")
    file.write("    current_model = target_model\n")
    file.write("    cmd.hide('all')\n")
    file.write("    cmd.color('grey60', 'elem C')\n")
    file.write("    cmd.show('sticks', f'resi {candidates[num-1][2]} and chain {candidates[num-1][3]}')\n")
    file.write("    cmd.color('cyan', 'elem C and visible')\n")
    file.write("    cmd.show('sticks', 'byres all within " + str(snakemake.config["search_dist"]) +
               " of (visible and sidechain)')\n")
    file.write("    label_b_factor(f'resi {candidates[num-1][2]} and chain {candidates[num-1][3]}', "
               "candidates[num-1][4])\n")
    file.write("    cmd.orient('visible')\n")
    file.write("    print(f'Showing candidate number {num}')\n\n\n")
    file.write("# Show the next candidate in the list.\n")
    file.write("def n():\n")
    file.write("    global can_num\n")
    file.write("    if can_num == len(candidates):\n")
    file.write("        can_num = 1\n")
    file.write("    else:\n")
    file.write("        can_num += 1\n")
    file.write("    goto(can_num)\n\n\n")
    file.write("# Show the previous candidate in the list.\n")
    file.write("def p():\n")
    file.write("    global can_num\n")
    file.write("    if can_num == 1:\n")
    file.write("        can_num = len(candidates)\n")
    file.write("    else:\n")
    file.write("        can_num -= 1\n")
    file.write("    goto(can_num)\n\n\n")
    file.write("# Print the number of candidates.\n")
    file.write("get_num()\n\n")
    file.write("# Load the first candidate.\n")
    file.write("goto(can_num)\n\n")
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
