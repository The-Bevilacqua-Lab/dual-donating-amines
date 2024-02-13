"""
This module will eventually do something.
"""

import sys
import os
import csv
import pandas as pd
from datetime import datetime
import const

os.chdir("/Users/drew/Documents/Research/Dual H-Bonding/Scratch/")
pd.set_option("display.width", 1000)
pd.set_option("display.max_columns", 40)
pd.set_option("display.max_rows", 1300)

# store a list of arguments provided when running this script
if len(sys.argv) > 1:
    args = sys.argv[1:]
else:
    args = []

# set the H-bonding criteria
# H_DIST_MAX = 3.5
# H_ANG_TOL = 180.0
H_DIST_MAX = 2.0
H_ANG_TOL = 30
DON_DIST_MAX = 3.5
DON_ANG_TOL = 45.0

# construct a list of tuples containing atom and residue names that describe atoms capable of donating and accepting an
# H-bond
don_acc_atoms = []
for residue in const.RESIDUE_LIBRARY:
    for donor in residue['don']:
        for acceptor in residue['acc']:
            if donor[0] == acceptor[0]:
                don_acc_atoms.append((residue['res'], donor[0]))

# extract the data from the hbond csv file and remove redundant lines
hbond_file = "NR_3.0_65061.22_hbond_test_data.csv"
# hbond_file = "NR_3.0_08591.1_hbond_test_data.csv"
hbond_col_names = ["don index", "don name", "don resn", "don resi", "don chain", "acc index", "acc name",
                   "acc resn", "acc resi", "acc chain", "dist", "ang", "vertex", "hydrogen",
                   "rotated side chain"]
# hbond_col_names = ["success", "don index", "don name", "don resn", "don resi", "don chain", "acc index", "acc name",
#                    "acc resn", "acc resi", "acc chain", "hbond", "dist", "ang", "vertex", "hydrogen",
#                    "rotated side chain"]
unique_col_comb = ["don index", "acc index", "hydrogen", "rotated side chain"]
hbond_data = pd.read_csv(hbond_file, names=hbond_col_names).drop_duplicates(subset=unique_col_comb)

# extract the data from the nuc csv file
nuc_file = "NR_3.0_08591.1_nuc_test_data.csv"
nuc_col_names = ["residue", "don index", "don name", "acc1 index", "acc1 name", "acc2 index", "acc2 name"]
nuc_data = pd.read_csv(nuc_file, names=nuc_col_names, index_col="residue")

# identify atom pairs that meet the H-bond criteria and include a donor of interest
don_hbonds = (hbond_data[(hbond_data["dist"] <= H_DIST_MAX) & (hbond_data["ang"] >= 180.0 - H_ANG_TOL)]
              .merge(pd.DataFrame(const.DONORS_OF_INTEREST, columns=["don resn", "don name"]), how='inner'))

# identify atom pairs that meet the H-bond criteria and include a protonated donor of interest
# values in the _merge column matching left_only indicate acceptors that cannot typically also donate an H-bond
prot_don_hbonds = (hbond_data[(hbond_data["dist"] <= H_DIST_MAX) & (hbond_data["ang"] >= 180.0 - H_ANG_TOL)]
                   .merge(pd.DataFrame(const.PROT_DONORS_OF_INTEREST, columns=["don resn", "don name"]), how='inner')
                   .merge(pd.DataFrame(don_acc_atoms, columns=["acc resn", "acc name"]), how='left', indicator=True))

# filter prot_don_hbonds such that only acceptors that cannot typically donate an H-bond are included
prot_don_hbonds_filtered = prot_don_hbonds[prot_don_hbonds["_merge"] == "left_only"].drop(columns="_merge")

# identify atom pairs that meet the H-bond criteria and include an acceptor of interest
acc_hbonds = (hbond_data[((hbond_data["vertex"] == "hydrogen") & (hbond_data["dist"] <= H_DIST_MAX) &
                          (hbond_data["ang"] >= 180.0 - H_ANG_TOL)) |
                         ((hbond_data["vertex"] == "donor") & (hbond_data["dist"] <= DON_DIST_MAX) &
                                    (hbond_data["ang"] >= 109.5 - DON_ANG_TOL) &
                                    (hbond_data["ang"] <= 109.5 + DON_ANG_TOL))]
              .merge(pd.DataFrame(const.ACCEPTORS_OF_INTEREST, columns=["acc resn", "acc name"]), how='inner'))

# Create a new list of H-bonding atom pairs involving the donors of interest with redundancy removed such that if
# multiple atom pairs involve the same donor and an acceptor belonging to different side chain conformations of the same
# ASN, GLN, or HIS residue, only the atom pairs with the original conformation are kept. There could be situations where
# one hydrogen of an exocyclic amine donates to the original side chain conformation of the acceptor residue while the
# other hydrogen of the same exocyclic amine donates to the rotated side chain conformation of the acceptor residue. For
# simplicity, H-bonds involving the original side chain conformations of acceptor residues are prioritized over
# situations where a dual H-bonding exocyclic amine could be inferred by discarding the atom pair involving the original
# side chain conformation and keeping the atom pair involving the rotated side chain conformation.
acc_grp = ["don index", "acc resn", "acc resi", "acc chain"]
don_hbonds_nonredundant = (don_hbonds[don_hbonds.groupby(acc_grp)["rotated side chain"]
                           .transform(lambda grp: grp.str.fullmatch("none")
                                      if any(grp.str.fullmatch("none")) else True)])

# Create a new list of H-bonding atom pairs involving the acceptors of interest with redundancy removed such that if
# multiple atom pairs involve the same acceptor and a donor belonging to different side chain conformations of the same
# ASN, GLN, or HIS residue, only the atom pairs with the original conformation are kept. Redundancy pertaining to
# donation by TYR OH atoms is also removed. The TYR OH hydrogen can occupy two different conformations. If both
# conformations are capable of forming an H-bond with the acceptor of interest, only the conformation with the shortest
# hydrogen-acceptor distance is kept.
don_grp = ["don resn", "don resi", "don chain", "acc index"]
acc_hbonds_nonredundant = (acc_hbonds[(((acc_hbonds["don resn"] != "TYR") &
                                        (acc_hbonds.groupby(don_grp)["rotated side chain"]
                                         .transform(lambda grp: grp.str.fullmatch("none")
                                                    if any(grp.str.fullmatch("none")) else True))) |
                                       ((acc_hbonds["don resn"] == "TYR") &
                                        (acc_hbonds.groupby(don_grp)["dist"]
                                         .transform(lambda grp: [mem == grp.min() for mem in grp]))))])

# identify nucleobases that donate either one or at least two H-bonds via their exocyclic amines
# the H-bonds from the latter nucleobases must involve both exocyclic amine hydrogens and at least two different
# acceptors
don_indices = list(don_hbonds_nonredundant.groupby("don index").groups.keys())
single_don_exo_amines = pd.Series(don_indices, index=don_indices, name="don index")[
    (don_hbonds_nonredundant.groupby("don index")["hydrogen"].nunique() == 1) |
    (don_hbonds_nonredundant.groupby("don index")["acc index"].nunique() == 1)]
dual_don_exo_amines = pd.Series(don_indices, index=don_indices, name="don index")[
    (don_hbonds_nonredundant.groupby("don index")["hydrogen"].nunique() == 2) &
    (don_hbonds_nonredundant.groupby("don index")["acc index"].nunique() >= 2)]

# if don_geom is provided as an argument, write H-bond geometry info to a csv file excluding certain atom pairs
if "don_geom" in args:
    don_acc_grp = ["don index", "acc index"]
    (don_hbonds_nonredundant[
        # do not include hydrogens with smaller D-H...A angle
        (don_hbonds_nonredundant.groupby(don_acc_grp)["ang"]
         .transform(lambda grp: [mem == grp.max() for mem in grp])) &
        # do not consider ASN, GLN, or HIS acceptors
        (~don_hbonds_nonredundant["acc resn"].isin(["ASN", "GLN", "HIS"])) &
        # do not consider H-bonds found in canonical WCF base pairs
        (~don_hbonds_nonredundant[["don name", "acc resn", "acc name"]].eq(["N6", "U", "O4"]).all(axis='columns')) &
        (~don_hbonds_nonredundant[["don name", "acc resn", "acc name"]].eq(["N4", "G", "O6"]).all(axis='columns')) &
        (~don_hbonds_nonredundant[["don name", "acc resn", "acc name"]].eq(["N2", "C", "O2"]).all(axis='columns'))]
     .to_csv("don_hbonds.csv", index=False, columns=["dist", "ang"]))

# if single_don_rev is provided as an argument, create a new python script with PyMOL commands to review the donors
# included in the single_don_exo_amines dataframe
if "single_don_rev" in args:
    with open(f"single_don_rev_{datetime.now().strftime('%y%m%d_%H%M%S_%f')}.py", "w") as file:
        file.write(f"# info on interactions came from {hbond_file}\n\n")
        file.write("from pymol import cmd\n\n")
        file.write("atoms = [\n")
        atom_dict = don_hbonds_nonredundant[
            # do not consider H-bonds found in canonical WCF base pairs
            (~don_hbonds_nonredundant[["don name", "acc resn", "acc name"]].eq(["N6", "U", "O4"]).all(axis='columns')) &
            (~don_hbonds_nonredundant[["don name", "acc resn", "acc name"]].eq(["N4", "G", "O6"]).all(axis='columns')) &
            (~don_hbonds_nonredundant[["don name", "acc resn", "acc name"]].eq(["N2", "C", "O2"]).all(axis='columns'))
        ].merge(single_don_exo_amines)[["hydrogen", "don resn", "don resi", "don chain", "acc name", "acc resn",
                                        "acc resi", "acc chain"]].to_dict('index')
        last_i = len(atom_dict) - 1
        for i, row in enumerate(atom_dict.values()):
            if i < last_i:
                file.write(f'    ["{row["hydrogen"]}", "{row["don resn"]}", "{row["don resi"]}", "{row["don chain"]}", '
                           f'"{row["acc name"]}", "{row["acc resn"]}", "{row["acc resi"]}", "{row["acc chain"]}"],\n')
            else:
                file.write(f'    ["{row["hydrogen"]}", "{row["don resn"]}", "{row["don resi"]}", "{row["don chain"]}", '
                           f'"{row["acc name"]}", "{row["acc resn"]}", "{row["acc resi"]}", "{row["acc chain"]}"]\n')
        file.write("]\n")
        file.write("\n")
        file.write("print(f'Number of atoms to review: {len(atoms)}')\n")
        file.write("\n")
        file.write("def num_atoms():\n")
        file.write("    print(len(atoms))\n")
        file.write("\n")
        file.write("def goto(num):\n")
        file.write("    i = int(num) - 1\n")
        file.write("    cmd.hide('all')\n")
        file.write("    cmd.color('grey60', 'elem C')\n")
        file.write("    cmd.show('sticks', f'resn {atoms[i][1]} and resi {atoms[i][2]} and chain {atoms[i][3]}')\n")
        file.write("    cmd.color('cyan', 'elem C and visible')\n")
        file.write("    cmd.show('sticks', 'byres all within 3.6 of (visible and sidechain)')\n")
        for don in const.PROT_DONORS_OF_INTEREST:
            file.write(f"    cmd.hide('sticks', '(neighbor (resn {don[0]} and name {don[1]})) and elem H')\n")
        file.write("    cmd.distance(f'dist_{num}', f'name {atoms[i][0]} and resn {atoms[i][1]} and resi {atoms[i][2]} "
                   "and chain {atoms[i][3]}', f'name {atoms[i][4]} and resn {atoms[i][5]} and resi {atoms[i][6]} and "
                   "chain {atoms[i][7]}')\n")
        file.write("    cmd.hide('labels')\n")
        file.write("    cmd.show('dashes', f'dist_{num}')\n")
        file.write("    cmd.orient('visible')\n")
        file.write("\n")
        file.write("cmd.extend('num_atoms', num_atoms)\n")
        file.write("cmd.extend('goto', goto)\n")
        file.write("\n")

# TODO prot donor should not donate to the same atom the exo amine is donating to (e.g., C788 A5 in 6xu8)

# inspect a specific region in the H-bond geometry
# print(don_hbonds_nonredundant[(don_hbonds_nonredundant["ang"] < 133) & (don_hbonds_nonredundant["ang"] > 132) & (don_hbonds_nonredundant["dist"] > 2.8) & (don_hbonds_nonredundant["dist"] < 2.85)][["dist", "ang"]])

# # print dual H-bonding nucleobases that also accept H-bonds via their acc of interest
# print(don_hbonds.merge(dual_don_exo_amines).merge(acc_hbonds_nonredundant, left_on=["don resn", "don resi", "don chain"], right_on=["acc resn", "acc resi", "acc chain"]))

# with open(f"acceptor_of_int_hbond_nonrot.csv", "w") as csv_file:
#     writer = csv.writer(csv_file)
#     for acceptor in h_bond_to_acc_nonredundant.keys():
#         for atom_pair_list in h_bond_to_acc_nonredundant[acceptor].values():
#             for atom_pair in atom_pair_list:
#                 if atom_pair["vertex"] == "hydrogen":
#                     writer.writerow([atom_pair["ang"], atom_pair["dist"]])
#
# with open(f"acceptor_of_int_hbond_rot.csv", "w") as csv_file:
#     writer = csv.writer(csv_file)
#     for acceptor in h_bond_to_acc_nonredundant.keys():
#         for atom_pair_list in h_bond_to_acc_nonredundant[acceptor].values():
#             for atom_pair in atom_pair_list:
#                 if atom_pair["vertex"] == "donor":
#                     writer.writerow([atom_pair["ang"], atom_pair["dist"]])

# count = 0
# for key in h_bond_to_acc_of_int.keys():
#     for x in h_bond_to_acc_of_int[key].values():
#         for y in x:
#             count += 1
# print(count)
