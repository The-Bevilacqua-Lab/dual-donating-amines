"""
Add a description for this module.
"""

from datetime import datetime
import residue_library


# Create a new python script with PyMOL commands to review the identified H-bonds. The function takes four arguments
# that specifies 1) the category of hbonds to consider, 2) the relevant pandas dataframe, 3) the name of the source
# file, and 4) the name of the output file.
def create_script(cat, df, source_files, output_file):
    with open(output_file, "w") as file:
        file.write("# this script is to be used with PyMOL\n")
        file.write("# info on interactions came from ")
        last_i = len(source_files) - 1
        for i, file_name in enumerate(source_files):
            if i == last_i:
                file.write("and ")
                file.write(file_name)
                file.write("\n\n")
            else:
                file.write(file_name)
                file.write(", ")
        file.write("from pymol import cmd\n\n\n")
        file.write("cmd.set('scene_animation_duration', '0.5')\n\n")
        if cat == "no_don_review":
            file.write("residues = [\n")
            for i in df.index.values:
                file.write(f'    ["{df.loc[i, "resn"]}", "{df.loc[i, "resi"]}", "{df.loc[i, "chain"]}"]')
                if i == df.index.values[-1]:
                    file.write('\n')
                else:
                    file.write(',\n')
            file.write("]\n")
        if cat == "single_don_review" or cat == "dual_don_review":
            file.write("hbonds = [\n")
            last_grp = list(df.groupby("don_index").groups.keys())[-1]
            for grp, mem in df.groupby("don_index"):
                for i, pair in enumerate(mem.to_dict('index').values()):
                    if i == 0:
                        file.write(f'    [["{pair["don_name"]}", "{pair["don_resn"]}", "{pair["don_resi"]}", '
                                   f'"{pair["don_chain"]}", "{pair["acc_name"]}", "{pair["acc_resn"]}", '
                                   f'"{pair["acc_resi"]}", "{pair["acc_chain"]}"]')
                    else:
                        file.write(f', ["{pair["don_name"]}", "{pair["don_resn"]}", "{pair["don_resi"]}", '
                                   f'"{pair["don_chain"]}", "{pair["acc_name"]}", "{pair["acc_resn"]}", '
                                   f'"{pair["acc_resi"]}", "{pair["acc_chain"]}"]')
                    if i == mem["don_index"].size - 1:
                        if grp != last_grp:
                            file.write('],\n')
                        else:
                            file.write(']\n')
            file.write("]\n")
        # # TODO the following if block needs review
        # if cat == "prot_don_review":
        #     file.write("hbonds = [\n")
        #     last_grp = list(df[0].groupby("don_index").groups.keys())[-1]
        #     for grp, mem in df[0].groupby("don_index"):
        #         for i, pair in enumerate(mem.to_dict('index').values()):
        #             if i == 0:
        #                 file.write(f'    [["{pair["hydrogen"]}", "{pair["don_resn"]}", "{pair["don_resi"]}", '
        #                            f'"{pair["don_chain"]}", "{pair["acc_name"]}", "{pair["acc_resn"]}", '
        #                            f'"{pair["acc_resi"]}", "{pair["acc_chain"]}"]')
        #             else:
        #                 file.write(f', ["{pair["hydrogen"]}", "{pair["don_resn"]}", "{pair["don_resi"]}", '
        #                            f'"{pair["don_chain"]}", "{pair["acc_name"]}", "{pair["acc_resn"]}", '
        #                            f'"{pair["acc_resi"]}", "{pair["acc_chain"]}"]')
        #             if i == mem["don_index"].size - 1:
        #                 if grp != last_grp:
        #                     file.write(f'],\n')
        #                 else:
        #                     file.write(f']\n')
        #     file.write("]\n")
        file.write("\n")
        if cat == "no_don_review":
            file.write("print(f'Number of residues to review: {len(residues)}')\n")
            file.write("\n")
            file.write("def num_residues():\n")
            file.write("    print(len(residues))\n")
        if cat == "single_don_review" or cat == "dual_don_review" or cat == "prot_don_review":
            file.write("print(f'Number of hbond sets to review: {len(hbonds)}')\n")
            file.write("\n")
            file.write("def num_hbonds():\n")
            file.write("    print(len(hbonds))\n")
        file.write("\n")
        if cat == "no_don_review":
            file.write("def goto(num):\n")
            file.write("    i = int(num) - 1\n")
            file.write("    cmd.hide('all')\n")
            file.write("    cmd.color('grey60', 'elem C')\n")
            file.write("    cmd.show('sticks', f'resn {residues[i][0]} and resi {residues[i][1]} and chain "
                       "{residues[i][2]}')\n")
            file.write("    cmd.color('cyan', 'elem C and visible')\n")
            file.write("    cmd.show('sticks', 'byres all within 4.1 of (visible and sidechain)')\n")
            for don in const.PROT_DONORS_OF_INTEREST:
                file.write(f"    cmd.hide('sticks', '(neighbor (resn {don[0]} and name {don[1]})) and elem H')\n")
            file.write("    cmd.orient('visible')\n")
        if cat == "single_don_review" or cat == "dual_don_review" or cat == "prot_don_review":
            file.write("def goto(num):\n")
            file.write("    i1 = int(num) - 1\n")
            file.write("    cmd.hide('all')\n")
            file.write("    cmd.color('grey60', 'elem C')\n")
            file.write("    cmd.show('sticks', f'resn {hbonds[i1][0][1]} and resi {hbonds[i1][0][2]} and chain "
                       "{hbonds[i1][0][3]}')\n")
            file.write("    cmd.color('cyan', 'elem C and visible')\n")
            file.write("    cmd.show('sticks', 'byres all within 4.1 of (visible and sidechain)')\n")
            for don in const.PROT_DONORS_OF_INTEREST:
                file.write(f"    cmd.hide('sticks', '(neighbor (resn {don[0]} and name {don[1]})) and elem H')\n")
            if cat == "prot_don_review":
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
        if cat == "no_don_review":
            file.write("for i in range(len(residues)):\n")
            file.write("    goto(i + 1)\n")
            file.write("    cmd.scene(f'{residues[i][0]}{residues[i][1]}_{residues[i][2]}', 'store')\n")
            file.write("\n")
            file.write("cmd.extend('num_residues', num_residues)\n")
        if cat == "single_don_review" or cat == "dual_don_review" or cat == "prot_don_review":
            file.write("for i in range(len(hbonds)):\n")
            file.write("    goto(i + 1)\n")
            file.write("    cmd.scene(f'{hbonds[i][0][1]}{hbonds[i][0][2]}_{hbonds[i][0][3]}', 'store')\n")
            file.write("\n")
            file.write("cmd.extend('num_hbonds', num_hbonds)\n")
        file.write("cmd.extend('goto', goto)\n")
        file.write("\n")
