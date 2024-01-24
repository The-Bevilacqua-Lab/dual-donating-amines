"""
This module will eventually do something.
"""

import sys
import os
import csv

os.chdir("/Users/drew/Documents/Research/Dual H-Bonding/Scratch")

# set the H-bonding criteria
# H_DIST_MAX = 2.5
# H_ANG_TOL = 60.0
H_DIST_MAX = 2.0
H_ANG_TOL = 30.0
DON_DIST_MAX = 3.5
DON_ANG_TOL = 45.0

# construct three tuples of tuples containing atom and residue names that describe either donor, protonated donor
# (formal charge of +1), or acceptor atoms of particular interest
DONORS_OF_INTEREST = (('A', 'N6'), ('C', 'N4'), ('G', 'N2'))
PROT_DONORS_OF_INTEREST = (('A', 'N1'), ('A', 'N3'), ('C', 'N3'), ('G', 'N3'))
ACCEPTORS_OF_INTEREST = (('A', 'N1'), ('A', 'N3'), ('C', 'O2'), ('C', 'N3'), ('G', 'N3'), ('G', 'O6'))

# extract the data from the csv file
with open(f"test_NR_3.0_08591.1_test_data.csv", "r") as csv_file:
    reader = csv.reader(csv_file)
    data = []
    for line in reader:
        data.append(
            {
                "don index": line[1],
                "don name": line[2],
                "don resn": line[3],
                "don resi": line[4],
                "don chain": line[5],
                "acc index": line[6],
                "acc name": line[7],
                "acc resn": line[8],
                "acc resi": line[9],
                "acc chain": line[10],
                "dist": line[12],
                "ang": line[13],
                "vertex": line[14].split(" ")[0],
                "hydrogen": line[15],
                "rotated side chain": line[16]
            })

# TODO remove redundancy from the different H's going to the same acceptor and from different conformation of the same acceptor
# identify nucleobases that donate H-bonds via their exocyclic amines
exocyclic_amine_donors = {}
for atom_pair in data:
    if (atom_pair["don resn"], atom_pair["don name"]) in DONORS_OF_INTEREST:
        if float(atom_pair["dist"]) <= H_DIST_MAX and float(atom_pair["ang"]) >= 180.0 - H_ANG_TOL:
            if atom_pair["don index"] not in exocyclic_amine_donors.keys():
                exocyclic_amine_donors[atom_pair["don index"]] = [atom_pair]
            else:
                exocyclic_amine_donors[atom_pair["don index"]].append(atom_pair)

# identify nucleobases that donate two H-bonds via their exocyclic amines
# each H-bond must involve a different hydrogen from the exocyclic amine and a different acceptor
# if an H-bond is donated to two different conformations of the same residue side chain, it only counts as one H-bond
dual_donation = {}
for exo_amine in exocyclic_amine_donors.values():
    if len(exo_amine) >= 2:
        filtered_pairs = []
        ambiguous_res = {}
        # add any H-bonding atom pairs that involve acceptors that belong to the side chains of ASN, GLN, and HIS
        # residues to a separate dictionary for evaluation in the subsequent for loop and add all other H-bonding atom
        # pairs to the filtered_pairs list
        for atom_pair in exo_amine:
            if atom_pair["acc resn"] in ["ASN", "GLN", "HIS"]:
                if atom_pair["acc resn"] + atom_pair["acc resi"] + atom_pair["acc chain"] not in ambiguous_res.keys():
                    ambiguous_res[atom_pair["acc resn"] + atom_pair["acc resi"] + atom_pair["acc chain"]] = [atom_pair]
                else:
                    ambiguous_res[atom_pair["acc resn"] + atom_pair["acc resi"] + atom_pair["acc chain"]].append(atom_pair)
            else:
                filtered_pairs.append(atom_pair)
        # if there are multiple H-bonding atom pairs where the acceptors belong to different side chain conformations of
        # the same ASN, GLN, or HIS residue, only add the H-bonding atom pairs with the acceptor that belongs to the
        # original conformation to the filtered_pairs list
        for res in ambiguous_res.values():
            original = []
            rotated = []
            for atom_pair in res:
                if atom_pair["rotated side chain"].split()[0] == "none":
                    original.append(atom_pair)
                elif atom_pair["rotated side chain"].split()[0] == "rotated":
                    rotated.append(atom_pair)
                else:
                    print(f"Error: The code was unable to categorize an acceptor in "
                          f"{atom_pair['acc resn']}{atom_pair['acc resi']} {atom_pair['acc chain']} as either original or "
                          f"rotated.")
                    sys.exit(1)
            if original:
                filtered_pairs.extend(original)
            elif rotated:
                filtered_pairs.extend(rotated)
        # TODO continue working with the code here
        h_name = []
        acc_indices = []
        for atom_pair in exo_amine:
            if atom_pair["hydrogen"] not in h_name:
                h_name.append(atom_pair["hydrogen"])
            if atom_pair["acc index"] not in acc_indices:
                acc_indices.append(atom_pair["acc index"])
        if len(h_name) == 2 and len(acc_indices) >= 2:
            print(exo_amine)

# count = 0
# for amine in exocyclic_amine_donors.values():
#     if len(amine) == 2:
#         print(amine[0]["acc name"])
#         print(amine[1]["acc name"])
#         print(amine[0]["don index"])
#         print()
#         count += 1
#
# print(count)
