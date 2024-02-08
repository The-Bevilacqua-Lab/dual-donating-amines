"""
This module will eventually do something.
"""

import sys
import os
import csv
import pandas as pd
import residue_library

pd.set_option("display.width", 500)
pd.set_option("display.max_columns", 17)
pd.set_option("display.max_rows", 1300)

os.chdir("/Users/drew/Documents/Research/Dual H-Bonding/Scratch/")

# set the H-bonding criteria
H_DIST_MAX = 3.6
H_ANG_TOL = 180.0
# H_DIST_MAX = 2.0
# H_ANG_TOL = 30
DON_DIST_MAX = 3.5
DON_ANG_TOL = 45.0

# construct three tuples of tuples containing atom and residue names that describe either donor, protonated donor
# (formal charge of +1), or acceptor atoms of particular interest
DONORS_OF_INTEREST = (('A', 'N6'), ('C', 'N4'), ('G', 'N2'))
PROT_DONORS_OF_INTEREST = (('A', 'N1'), ('A', 'N3'), ('C', 'N3'), ('G', 'N3'))
ACCEPTORS_OF_INTEREST = (('A', 'N1'), ('A', 'N3'), ('C', 'O2'), ('C', 'N3'), ('G', 'N3'), ('G', 'O6'))

# construct a list of tuples containing atom and residue names that describe atoms capable of donating and accepting an
# H-bond
don_acc_atoms = []
for residue in residue_library.residue_library:
    for donor in residue['don']:
        for acceptor in residue['acc']:
            if donor[0] == acceptor[0]:
                don_acc_atoms.append((residue['res'], donor[0]))

# extract the data from the hbond csv file and remove redundant lines
hbond_file = "NR_3.0_08591.1_hbond_test_data.csv"
hbond_col_names = ["success", "don index", "don name", "don resn", "don resi", "don chain", "acc index", "acc name",
                   "acc resn", "acc resi", "acc chain", "hbond", "dist", "ang", "vertex", "hydrogen",
                   "rotated side chain"]
unique_col_comb = ["don index", "acc index", "hydrogen", "rotated side chain"]
hbond_data = pd.read_csv(hbond_file, names=hbond_col_names).drop_duplicates(subset=unique_col_comb)

# extract the data from the nuc csv file
nuc_file = "NR_3.0_08591.1_nuc_test_data.csv"
nuc_col_names = ["residue", "don index", "don name", "acc1 index", "acc1 name", "acc2 index", "acc2 name"]
nuc_data = pd.read_csv(nuc_file, names=nuc_col_names, index_col="residue")

# identify atom pairs that meet the H-bond criteria and include a donor of interest
don_hbonds = (hbond_data[(hbond_data["dist"] <= H_DIST_MAX) & (hbond_data["ang"] >= 180.0 - H_ANG_TOL)]
              .merge(pd.DataFrame(DONORS_OF_INTEREST, columns=["don resn", "don name"]), how='inner'))

# identify atom pairs that meet the H-bond criteria and include a protonated donor of interest
# values in the _merge column matching left_only indicate acceptors that cannot typically also donate an H-bond
prot_don_hbonds = (hbond_data[(hbond_data["dist"] <= H_DIST_MAX) & (hbond_data["ang"] >= 180.0 - H_ANG_TOL)]
                   .merge(pd.DataFrame(PROT_DONORS_OF_INTEREST, columns=["don resn", "don name"]), how='inner')
                   .merge(pd.DataFrame(don_acc_atoms, columns=["acc resn", "acc name"]), how='left', indicator=True))

# filter prot_don_hbonds such that only acceptors that cannot typically donate an H-bond are included
prot_don_hbonds_filtered = prot_don_hbonds[prot_don_hbonds["_merge"] == "left_only"]

# identify atom pairs that meet the H-bond criteria and include an acceptor of interest
acc_hbonds = (hbond_data[((hbond_data["vertex"] == "hydrogen") & (hbond_data["dist"] <= H_DIST_MAX) &
                          (hbond_data["ang"] >= 180.0 - H_ANG_TOL)) |
                         ((hbond_data["vertex"] == "donor") & (hbond_data["dist"] <= DON_DIST_MAX) &
                                    (hbond_data["ang"] >= 109.5 - DON_ANG_TOL) &
                                    (hbond_data["ang"] <= 109.5 + DON_ANG_TOL))]
              .merge(pd.DataFrame(ACCEPTORS_OF_INTEREST, columns=["acc resn", "acc name"]), how='inner'))

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
single_don_exo_amines = pd.Series(don_indices, index=don_indices)[
    (don_hbonds_nonredundant.groupby("don index")["hydrogen"].nunique() == 1) |
    (don_hbonds_nonredundant.groupby("don index")["acc index"].nunique() == 1)]
dual_don_exo_amines = pd.Series(don_indices, index=don_indices)[
    (don_hbonds_nonredundant.groupby("don index")["hydrogen"].nunique() == 2) &
    (don_hbonds_nonredundant.groupby("don index")["acc index"].nunique() >= 2)]

don_acc_grp = ["don index", "acc index"]
# print(don_hbonds_nonredundant[
#           # do not include hydrogens with smaller D-H...A angle
#           (don_hbonds_nonredundant.groupby(don_acc_grp)["ang"]
#            .transform(lambda grp: [mem == grp.max() for mem in grp])) &
#           # do not consider ASN, GLN, or HIS acceptors
#           (~don_hbonds_nonredundant["acc resn"].isin(["ASN", "GLN", "HIS"])) &
#           ((don_hbonds_nonredundant["don name"] != "N2") & (don_hbonds_nonredundant["acc resn", "acc name"] != ("C", "O2")))]
#       .to_csv("don_hbonds.csv", index=False, columns=["dist", "ang"]))

print(don_hbonds_nonredundant[((don_hbonds_nonredundant["don name"] != "N2") & (don_hbonds_nonredundant["acc resn"] != "C") & (don_hbonds_nonredundant["acc name"] != "O2"))])

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
