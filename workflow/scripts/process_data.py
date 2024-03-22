"""
This module will eventually do something.
"""

import sys
import os
import pandas as pd
import csv
import const
import review_hbonds

pd.set_option("display.width", 1000)
pd.set_option("display.max_columns", 40)
pd.set_option("display.max_rows", 1300)

# set the H-bonding criteria
# H_DIST_MAX = 20
# H_ANG_TOL = 1000
H_DIST_MAX = 2.0
H_ANG_TOL = 30
DON_DIST_MAX = 3.5
DON_ANG_TOL = 45.0

# construct a list of tuples containing atom and residue names that describe atoms capable of both donating and
# accepting an H-bond
don_acc_atoms = []
for residue in const.RESIDUE_LIBRARY:
    for donor in residue['don']:
        for acceptor in residue['acc']:
            if donor[0] == acceptor[0]:
                don_acc_atoms.append((residue['res'], donor[0]))

# extract the data from the hbond csv file and remove redundant lines
hbond_col_names = ["don index", "don name", "don resn", "don resi", "don chain", "acc index", "acc name",
                   "acc resn", "acc resi", "acc chain", "dist", "ang", "vertex", "hydrogen",
                   "rotated side chain"]
unique_col_comb = ["don index", "acc index", "hydrogen", "rotated side chain"]
hbond_data = (pd.read_csv(snakemake.input.hbond, names=hbond_col_names, comment="#")
              .drop_duplicates(subset=unique_col_comb))

# extract the data from the nuc csv file
nuc_col_names = ["index", "atom name", "resn", "resi", "chain"]
nuc_data = pd.read_csv(snakemake.input.nuc, names=nuc_col_names, comment="#")

# identify atom pairs that meet the H-bond criteria and include a donor of interest
don_hbonds = (hbond_data[(hbond_data["dist"] <= H_DIST_MAX) & (hbond_data["ang"] >= 180.0 - H_ANG_TOL)]
              .merge(pd.DataFrame(const.DONORS_OF_INTEREST, columns=["don resn", "don name"]), how='inner'))

# identify atom pairs that meet the H-bond criteria and include a protonated donor of interest
# values in the _merge column matching left_only indicate acceptors that cannot typically also donate an H-bond
prot_don_hbonds_unfiltered = (hbond_data[(hbond_data["dist"] <= H_DIST_MAX) & (hbond_data["ang"] >= 180.0 - H_ANG_TOL)]
                              .merge(pd.DataFrame(const.PROT_DONORS_OF_INTEREST, columns=["don resn", "don name"]),
                                     how='inner')
                              .merge(pd.DataFrame(don_acc_atoms, columns=["acc resn", "acc name"]), how='left',
                                     indicator=True))

# Filter prot_don_hbonds_unfiltered such that only acceptors that cannot typically also donate an H-bond are included.
# Also remove side chain acceptor atoms from ASN and GLN residues due to ambiguity in side chain conformation.
prot_don_hbonds = (prot_don_hbonds_unfiltered[(prot_don_hbonds_unfiltered["_merge"] == "left_only") &
                                              ~((prot_don_hbonds_unfiltered["acc resn"] == "ASN") &
                                                (prot_don_hbonds_unfiltered["acc name"] == "OD1")) &
                                              ~((prot_don_hbonds_unfiltered["acc resn"] == "GLN") &
                                                (prot_don_hbonds_unfiltered["acc name"] == "OE1"))]
                   .drop(columns="_merge"))

# identify atom pairs that meet the H-bond criteria and include an acceptor of interest
acc_hbonds = (hbond_data[((hbond_data["vertex"] == "hydrogen") & (hbond_data["dist"] <= H_DIST_MAX) &
                          (hbond_data["ang"] >= 180.0 - H_ANG_TOL)) |
                         ((hbond_data["vertex"] == "donor") & (hbond_data["dist"] <= DON_DIST_MAX) &
                                    (hbond_data["ang"] >= 109.5 - DON_ANG_TOL) &
                                    (hbond_data["ang"] <= 109.5 + DON_ANG_TOL))]
              .merge(pd.DataFrame(const.ACCEPTORS_OF_INTEREST, columns=["acc resn", "acc name"]), how='inner'))

# identify atom pairs that meet the H-bond criteria and include a deprotonated acceptor of interest
# values in the _merge column matching left_only indicate donors that cannot typically also accept an H-bond
deprot_acc_hbonds_unfiltered = (hbond_data[((hbond_data["vertex"] == "hydrogen") & (hbond_data["dist"] <= H_DIST_MAX) &
                                            (hbond_data["ang"] >= 180.0 - H_ANG_TOL)) |
                                           ((hbond_data["vertex"] == "donor") & (hbond_data["dist"] <= DON_DIST_MAX) &
                                            (hbond_data["ang"] >= 109.5 - DON_ANG_TOL) &
                                            (hbond_data["ang"] <= 109.5 + DON_ANG_TOL))]
                                .merge(pd.DataFrame(const.DEPROT_ACCEPTORS_OF_INTEREST,
                                       columns=["don resn", "don name"]), how='inner')
                                .merge(pd.DataFrame(don_acc_atoms, columns=["don resn", "don name"]), how='left',
                                       indicator=True))

# Filter deprot_acc_hbonds_unfiltered such that only donors that cannot typically also accept an H-bond are included.
# Also remove side chain donor atoms from ASN and GLN residues due to ambiguity in side chain conformation.
deprot_acc_hbonds = (deprot_acc_hbonds_unfiltered[(deprot_acc_hbonds_unfiltered["_merge"] == "left_only") &
                                                  ~((deprot_acc_hbonds_unfiltered["acc resn"] == "ASN") &
                                                    (deprot_acc_hbonds_unfiltered["acc name"] == "ND2")) &
                                                  ~((deprot_acc_hbonds_unfiltered["acc resn"] == "GLN") &
                                                    (deprot_acc_hbonds_unfiltered["acc name"] == "NE2"))]
                     .drop(columns="_merge"))

# Create a new list of H-bonding atom pairs involving the donors of interest with redundancy removed such that if
# multiple atom pairs involve the same donor and an acceptor belonging to different side chain conformations of the same
# ASN, GLN, or HIS residue, only the atom pairs with the original conformation are kept. There could be situations where
# one hydrogen of an exocyclic amine donates to the original side chain conformation of the acceptor residue while the
# other hydrogen of the same exocyclic amine donates to the rotated side chain conformation of the acceptor residue. For
# simplicity, H-bonds involving the original side chain conformations of acceptor residues are prioritized over
# situations where a dual H-bonding exocyclic amine could be inferred by discarding the atom pair involving the original
# side chain conformation and keeping the atom pair involving the rotated side chain conformation.
acc_grp = ["don index", "acc resn", "acc resi", "acc chain"]
don_hbonds_nr = (don_hbonds[don_hbonds.groupby(acc_grp)["rotated side chain"]
                 .transform(lambda grp: grp.str.fullmatch("none") if any(grp.str.fullmatch("none")) else True)])

# Create a new list of H-bonding atom pairs involving the acceptors of interest with redundancy removed such that if
# multiple atom pairs involve the same acceptor and a donor belonging to different side chain conformations of the same
# ASN, GLN, or HIS residue, only the atom pairs with the original conformation are kept. Redundancy pertaining to
# donation by TYR OH atoms is also removed. The TYR OH hydrogen can occupy two different conformations. If both
# conformations are capable of forming an H-bond with the acceptor of interest, only the conformation with the shortest
# hydrogen-acceptor distance is kept.
don_grp = ["don resn", "don resi", "don chain", "acc index"]
acc_hbonds_nr = (acc_hbonds[(((acc_hbonds["don resn"] != "TYR") &
                              (acc_hbonds.groupby(don_grp)["rotated side chain"]
                               .transform(lambda grp: grp.str.fullmatch("none")
                                          if any(grp.str.fullmatch("none")) else True))) |
                             ((acc_hbonds["don resn"] == "TYR") &
                              (acc_hbonds.groupby(don_grp)["dist"]
                               .transform(lambda grp: [mem == grp.min() for mem in grp]))))])

# identify nucleobases that donate either one or at least two H-bonds via their exocyclic amines
# the H-bonds from the latter nucleobases must involve both exocyclic amine hydrogens and at least two different
# acceptors
if don_hbonds_nr.size > 0:
    don_indices = list(don_hbonds_nr.groupby("don index").groups.keys())
    single_don_exo_amines = pd.Series(don_indices, index=don_indices, name="don index")[
        (don_hbonds_nr.groupby("don index")["hydrogen"].nunique() == 1) |
        (don_hbonds_nr.groupby("don index")["acc index"].nunique() == 1)]
    dual_don_exo_amines = pd.Series(don_indices, index=don_indices, name="don index")[
        (don_hbonds_nr.groupby("don index")["hydrogen"].nunique() == 2) &
        (don_hbonds_nr.groupby("don index")["acc index"].nunique() >= 2)]
else:
    single_don_exo_amines = pd.Series([])
    dual_don_exo_amines = pd.Series([])

# write dataframes to csv files
don_hbonds_nr.to_csv(snakemake.output.don_hbonds_nr)
prot_don_hbonds.to_csv(snakemake.output.prot_don_hbonds)
acc_hbonds_nr.to_csv(snakemake.output.acc_hbonds_nr)
deprot_acc_hbonds.to_csv(snakemake.output.deprot_acc_hbonds)
single_don_exo_amines.to_csv(snakemake.output.single_don_exo_amines)
dual_don_exo_amines.to_csv(snakemake.output.dual_don_exo_amines)

# write H-bond geometry info related to don_hbonds_nr excluding certain atom pairs to a csv file
if don_hbonds_nr.size > 0:
    don_acc_grp = ["don index", "acc index"]
    (don_hbonds_nr[
         # only include hydrogens with smaller D-H...A distances
         (don_hbonds_nr.groupby(don_acc_grp)["dist"]
          .transform(lambda grp: [mem == grp.min() for mem in grp])) &
         # do not consider H-bonds found in canonical WCF base pairs
         (~don_hbonds_nr[["don name", "acc resn", "acc name"]].eq(["N6", "U", "O4"]).all(axis='columns')) &
         (~don_hbonds_nr[["don name", "acc resn", "acc name"]].eq(["N4", "G", "O6"]).all(axis='columns')) &
         (~don_hbonds_nr[["don name", "acc resn", "acc name"]].eq(["N2", "C", "O2"]).all(axis='columns'))]
     .to_csv(snakemake.output.don_hbonds_geom, index=False, columns=["dist", "ang"]))
else:
    with open(snakemake.output.don_hbonds_geom, "w") as write_file:
        csv.writer(write_file).writerow(["dist", "ang"])

# # write H-bond geometry info related to prot_don_hbonds to a csv file
# prot_don_hbonds.to_csv(snakemake.output.prot_don_hbonds, index=False, columns=["dist", "ang"])
#
# # create a new python script with PyMOL commands to review the donors included in the single_don_exo_amines dataframe
# review_hbonds.create_script("single_don_rev", [don_hbonds_nr, single_don_exo_amines], snakemake.input.hbond, snakemake.output.single_don_rev)
#
# # create a new python script with PyMOL commands to review the donors included in the dual_don_exo_amines dataframe
# review_hbonds.create_script("dual_don_rev", [don_hbonds_nr, dual_don_exo_amines], snakemake.input.hbond, snakemake.output.dual_don_rev)
#
# # create a new python script with PyMOL commands to review the donors included in the prot_don_hbonds dataframe
# review_hbonds.create_script("prot_don_rev", [prot_don_hbonds], snakemake.input.hbond, snakemake.output.prot_don_rev)

# TODO prot donor should not donate to the same atom the exo amine is donating to (e.g., C788 A5 in 6xu8)
# TODO keep in mind ["H04", "A", "4212", "L5", "OD1", "ASN", "3", "LT"] in 8GLP
# TODO why is NR_3.0_65061.22_hbond_test_data/LA/LA/ARG`3/NH2 missing a proton?

# inspect a specific region in the H-bond geometry
# print(don_hbonds_nr[(don_hbonds_nr["ang"] < 133) & (don_hbonds_nr["ang"] > 132) & (don_hbonds_nr["dist"] > 2.8) & (don_hbonds_nr["dist"] < 2.85)][["dist", "ang"]])

# # print dual H-bonding nucleobases that also accept H-bonds via their acc of interest
# print(don_hbonds.merge(dual_don_exo_amines).merge(acc_hbonds_nr, left_on=["don resn", "don resi", "don chain"], right_on=["acc resn", "acc resi", "acc chain"]))

# count = 0
# for key in h_bond_to_acc_of_int.keys():
#     for x in h_bond_to_acc_of_int[key].values():
#         for y in x:
#             count += 1
# print(count)
