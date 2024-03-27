"""
This module will eventually do something.
"""

from datetime import datetime
import const


# Create a new python script with PyMOL commands to review the identified H-bonds. The function takes three arguments
# that 1) specifies the category of hbonds to consider, 2) provides a list of the relevant pandas dataframes, and
# 3) provides the name of the source file.
def create_script(cat, df_list, source_file, output_file):
    with open(output_file, "w") as file:
        file.write("# this script is to be used with PyMOL\n")
        file.write(f"# info on interactions came from {source_file}\n")
        file.write("# for best results, use with modified .cif file created by collect_data.py\n\n")
        file.write("from pymol import cmd\n\n\n")
        file.write("cmd.set('scene_animation_duration', '0.5')\n\n")
        if cat == "single_don_rev" or cat == "dual_don_rev":
            file.write("hbonds = [\n")
            last_grp = list(df_list[0].merge(df_list[1]).groupby("don_index").groups.keys())[-1]
            for grp, mem in df_list[0].merge(df_list[1]).groupby("don_index"):
                for i, pair in enumerate(mem.to_dict('index').values()):
                    if i == 0:
                        file.write(f'    [["{pair["hydrogen"]}", "{pair["don_resn"]}", "{pair["don_resi"]}", '
                                   f'"{pair["don_chain"]}", "{pair["acc_name"]}", "{pair["acc_resn"]}", '
                                   f'"{pair["acc_resi"]}", "{pair["acc_chain"]}"]')
                    else:
                        file.write(f', ["{pair["hydrogen"]}", "{pair["don_resn"]}", "{pair["don_resi"]}", '
                                   f'"{pair["don_chain"]}", "{pair["acc_name"]}", "{pair["acc_resn"]}", '
                                   f'"{pair["acc_resi"]}", "{pair["acc_chain"]}"]')
                    if i == mem["don_index"].size - 1:
                        if grp != last_grp:
                            file.write(f'],\n')
                        else:
                            file.write(f']\n')
            file.write("]\n")
        if cat == "prot_don_rev":
            file.write("hbonds = [\n")
            last_grp = list(df_list[0].groupby("don_index").groups.keys())[-1]
            for grp, mem in df_list[0].groupby("don_index"):
                for i, pair in enumerate(mem.to_dict('index').values()):
                    if i == 0:
                        file.write(f'    [["{pair["hydrogen"]}", "{pair["don_resn"]}", "{pair["don_resi"]}", '
                                   f'"{pair["don_chain"]}", "{pair["acc_name"]}", "{pair["acc_resn"]}", '
                                   f'"{pair["acc_resi"]}", "{pair["acc_chain"]}"]')
                    else:
                        file.write(f', ["{pair["hydrogen"]}", "{pair["don_resn"]}", "{pair["don_resi"]}", '
                                   f'"{pair["don_chain"]}", "{pair["acc_name"]}", "{pair["acc_resn"]}", '
                                   f'"{pair["acc_resi"]}", "{pair["acc_chain"]}"]')
                    if i == mem["don_index"].size - 1:
                        if grp != last_grp:
                            file.write(f'],\n')
                        else:
                            file.write(f']\n')
            file.write("]\n")
        file.write("\n")
        file.write("print(f'Number of hbond sets to review: {len(hbonds)}')\n")
        file.write("\n")
        file.write("def num_hbonds():\n")
        file.write("    print(len(hbonds))\n")
        file.write("\n")
        if cat == "single_don_rev" or cat == "dual_don_rev" or cat == "prot_don_rev":
            file.write("def goto(num):\n")
            file.write("    i1 = int(num) - 1\n")
            file.write("    cmd.hide('all')\n")
            file.write("    cmd.color('grey60', 'elem C')\n")
            file.write("    cmd.show('sticks', f'resn {hbonds[i1][0][1]} and resi {hbonds[i1][0][2]} and chain "
                       "{hbonds[i1][0][3]}')\n")
            file.write("    cmd.color('cyan', 'elem C and visible')\n")
            file.write("    cmd.show('sticks', 'byres all within 3.6 of (visible and sidechain)')\n")
            for don in const.PROT_DONORS_OF_INTEREST:
                file.write(f"    cmd.hide('sticks', '(neighbor (resn {don[0]} and name {don[1]})) and elem H')\n")
            if cat == "prot_don_rev":
                file.write("    cmd.show('sticks', f'name {hbonds[i1][0][0]} and resn {hbonds[i1][0][1]} and "
                           "resi {hbonds[i1][0][2]} and chain {hbonds[i1][0][3]}')\n")
            file.write("    for i2, pair in enumerate(hbonds[i1]):\n")
            file.write("        cmd.distance(f'dist_{num}_{i2}', f'name {hbonds[i1][i2][0]} and "
                       "resn {hbonds[i1][i2][1]} and resi {hbonds[i1][i2][2]} and chain {hbonds[i1][i2][3]}', "
                       "f'name {hbonds[i1][i2][4]} and resn {hbonds[i1][i2][5]} and resi {hbonds[i1][i2][6]} and "
                       "chain {hbonds[i1][i2][7]}')\n")
            file.write("    cmd.hide('labels')\n")
            file.write("    cmd.show('dashes', f'dist_{num}*')\n")
            file.write("    cmd.orient('visible')\n")
        file.write("\n")
        file.write("for i in range(len(hbonds)):\n")
        file.write("    goto(i)\n")
        file.write("    cmd.scene(str(i), 'store')\n")
        file.write("\n")
        file.write("cmd.extend('num_hbonds', num_hbonds)\n")
        file.write("cmd.extend('goto', goto)\n")
        file.write("\n")
