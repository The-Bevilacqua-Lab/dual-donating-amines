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
H_DIST_MAX = 3.5
H_ANG_TOL = 180
DON_DIST_MAX = 3.5
DON_ANG_TOL = 45.0

# construct three tuples of tuples containing atom and residue names that describe either donor, protonated donor
# (formal charge of +1), or acceptor atoms of particular interest
DONORS_OF_INTEREST = (('A', 'N6'), ('C', 'N4'), ('G', 'N2'))
PROT_DONORS_OF_INTEREST = (('A', 'N1'), ('A', 'N3'), ('C', 'N3'), ('G', 'N3'))
ACCEPTORS_OF_INTEREST = (('A', 'N1'), ('A', 'N3'), ('C', 'O2'), ('C', 'N3'), ('G', 'N3'), ('G', 'O6'))

# extract the data from the hbond csv file
hbond_file = "NR_3.0_65061.22_hbond_test_data.csv"
with open(hbond_file, "r") as csv_file:
    reader = csv.reader(csv_file)
    hbond_data = []
    for line in reader:
        hbond_data.append(
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
                "vertex": line[14],
                "hydrogen": line[15],
                "rotated side chain": line[16]
            })

# extract the data from the nuc csv file
nuc_file = "NR_3.0_65061.22_nuc_test_data.csv"
with open(nuc_file, "r") as csv_file:
    reader = csv.reader(csv_file)
    nuc_ids = []
    nuc_data = {}
    for line in reader:
        if line[0] not in nuc_ids:
            nuc_ids.append(line[0])
            nuc_data[line[0]] = ((line[1], line[2]), (line[3], line[4]), (line[5], line[6]))
        else:
            print(f"Error: There are multiple lines referencing the same nucleobase in {nuc_file}.")
            sys.exit(1)

# identify atom pairs that meet the H-bond criteria and include a donor of interest
# also, identify atom pairs that meet the H-bond criteria and include an acceptor of interest
h_bond_from_don_of_int = {}
h_bond_to_acc_of_int = {
    "A(N1)": {},
    "A(N3)": {},
    "C(O2)": {},
    "C(N3)": {},
    "G(N3)": {},
    "G(O6)": {}
}
for atom_pair in hbond_data:
    # for exocyclic amines donors
    if (atom_pair["don resn"], atom_pair["don name"]) in DONORS_OF_INTEREST:
        if float(atom_pair["dist"]) <= H_DIST_MAX and float(atom_pair["ang"]) >= 180.0 - H_ANG_TOL:
            don_nucleobase = atom_pair["don resn"] + atom_pair["don resi"] + atom_pair["don chain"]
            if don_nucleobase not in h_bond_from_don_of_int.keys():
                h_bond_from_don_of_int[don_nucleobase] = [atom_pair]
            else:
                already_included = False
                for included_atom_pair in h_bond_from_don_of_int[don_nucleobase]:
                    if (atom_pair["don index"] == included_atom_pair["don index"] and
                            atom_pair["acc index"] == included_atom_pair["acc index"] and
                            atom_pair["hydrogen"] == included_atom_pair["hydrogen"] and
                            atom_pair["rotated side chain"] == included_atom_pair["rotated side chain"]):
                        already_included = True
                if not already_included:
                    h_bond_from_don_of_int[don_nucleobase].append(atom_pair)
    # for acceptors of interest
    if (atom_pair["acc resn"], atom_pair["acc name"]) in ACCEPTORS_OF_INTEREST:
        if ((atom_pair["vertex"] == "hydrogen" and float(atom_pair["dist"]) <= H_DIST_MAX and
                float(atom_pair["ang"]) >= 180.0 - H_ANG_TOL) or
                (atom_pair["vertex"] == "donor" and float(atom_pair["dist"]) <= DON_DIST_MAX and
                 109.5 - DON_ANG_TOL <= float(atom_pair["ang"]) <= 109.5 + DON_ANG_TOL)):
            acc_nucleobase = atom_pair["acc resn"] + atom_pair["acc resi"] + atom_pair["acc chain"]
            acc_of_int = f'{atom_pair["acc resn"]}({atom_pair["acc name"]})'
            if acc_nucleobase not in h_bond_to_acc_of_int[acc_of_int].keys():
                h_bond_to_acc_of_int[acc_of_int][acc_nucleobase] = [atom_pair]
            else:
                already_included = False
                for included_atom_pair in h_bond_to_acc_of_int[acc_of_int][acc_nucleobase]:
                    if (atom_pair["don index"] == included_atom_pair["don index"] and
                            atom_pair["acc index"] == included_atom_pair["acc index"] and
                            atom_pair["hydrogen"] == included_atom_pair["hydrogen"] and
                            atom_pair["rotated side chain"] == included_atom_pair["rotated side chain"]):
                        already_included = True
                if not already_included:
                    h_bond_to_acc_of_int[acc_of_int][acc_nucleobase].append(atom_pair)

# create a new list of H-bonding atom pairs involving the donors of interest with redundancy removed such that if
# multiple atom pairs involve different side chain conformations of the same ASN, GLN, or HIS residue, only the atom
# pairs with the original conformation are kept
h_bond_from_don_nonredundant = {}
for donor in h_bond_from_don_of_int.values():
    don_nucleobase = donor[0]["don resn"] + donor[0]["don resi"] + donor[0]["don chain"]
    filtered_pairs = []
    ambiguous_res = {}
    # add any H-bonding atom pairs that involve acceptors that belong to the side chains of ASN, GLN, and HIS
    # residues to a separate dictionary for evaluation in the subsequent for loop and add all other H-bonding atom
    # pairs to the filtered_pairs list
    for atom_pair in donor:
        if atom_pair["acc resn"] in ["ASN", "GLN", "HIS"]:
            acc_res = atom_pair["acc resn"] + atom_pair["acc resi"] + atom_pair["acc chain"]
            if acc_res not in ambiguous_res.keys():
                ambiguous_res[acc_res] = [atom_pair]
            else:
                ambiguous_res[acc_res].append(atom_pair)
        else:
            filtered_pairs.append(atom_pair)
    # if there are multiple H-bonding atom pairs where the acceptors belong to different side chain conformations of
    # the same ASN, GLN, or HIS residue, only add the H-bonding atom pairs with the acceptor that belongs to the
    # original conformation to the filtered_pairs list
    for res in ambiguous_res.values():
        original = []
        rotated = []
        for atom_pair in res:
            if atom_pair["rotated side chain"] == "none":
                original.append(atom_pair)
            else:
                rotated.append(atom_pair)
        if original:
            filtered_pairs.extend(original)
        elif rotated:
            filtered_pairs.extend(rotated)
    h_bond_from_don_nonredundant[don_nucleobase] = filtered_pairs

# create a new list of H-bonding atom pairs involving the acceptors of interest with redundancy removed such that if
# multiple atom pairs involve different side chain conformations of the same ASN, GLN, or HIS residue, only the atom
# pairs with the original conformation are kept
h_bond_to_acc_nonredundant = {
    "A(N1)": {},
    "A(N3)": {},
    "C(O2)": {},
    "C(N3)": {},
    "G(N3)": {},
    "G(O6)": {}
}
for acceptor in h_bond_to_acc_of_int.keys():
    for atom_pair_list in h_bond_to_acc_of_int[acceptor].values():
        acc_nucleobase = atom_pair_list[0]["acc resn"] + atom_pair_list[0]["acc resi"] + atom_pair_list[0]["acc chain"]
        filtered_pairs = []
        ambiguous_res = {}
        # add any H-bonding atom pairs that involve donors that belong to the side chains of ASN, GLN, and HIS
        # residues to a separate dictionary for evaluation in the subsequent for loop and add all other H-bonding atom
        # pairs to the filtered_pairs list
        for atom_pair in atom_pair_list:
            if atom_pair["don resn"] in ["ASN", "GLN", "HIS"]:
                don_res = atom_pair["don resn"] + atom_pair["don resi"] + atom_pair["don chain"]
                if don_res not in ambiguous_res.keys():
                    ambiguous_res[don_res] = [atom_pair]
                else:
                    ambiguous_res[don_res].append(atom_pair)
            else:
                filtered_pairs.append(atom_pair)
        # if there are multiple H-bonding atom pairs where the donors belong to different side chain conformations of
        # the same ASN, GLN, or HIS residue, only add the H-bonding atom pairs with the donor that belongs to the
        # original conformation to the filtered_pairs list
        for res in ambiguous_res.values():
            original = []
            rotated = []
            for atom_pair in res:
                if atom_pair["rotated side chain"] == "none":
                    original.append(atom_pair)
                else:
                    rotated.append(atom_pair)
            if original:
                filtered_pairs.extend(original)
            elif rotated:
                filtered_pairs.extend(rotated)
        h_bond_to_acc_nonredundant[acceptor][acc_nucleobase] = filtered_pairs

# identify nucleobases that donate either one or at least two H-bonds via their exocyclic amines
# the H-bonds from the latter nucleobases must involve both exocyclic amine hydrogens and at least two different
# acceptors
single_donating_exo_amines = {}
dual_donating_exo_amines = {}
for donor in h_bond_from_don_nonredundant.values():
    don_nucleobase = donor[0]["don resn"] + donor[0]["don resi"] + donor[0]["don chain"]
    if len(donor) == 1:
        single_donating_exo_amines[don_nucleobase] = donor
    elif len(donor) >= 2:
        h_name = []
        acc_indices = []
        for atom_pair in donor:
            if atom_pair["hydrogen"] not in h_name:
                h_name.append(atom_pair["hydrogen"])
            if atom_pair["acc index"] not in acc_indices:
                acc_indices.append(atom_pair["acc index"])
        if len(h_name) == 2 and len(acc_indices) >= 2:
            dual_donating_exo_amines[don_nucleobase] = donor
        else:
            single_donating_exo_amines[don_nucleobase] = donor


with open(f"donor_of_int_hbond.csv", "w") as csv_file:
    writer = csv.writer(csv_file)
    for atom_pair_list in h_bond_from_don_nonredundant.values():
        keep = []
        for atom_pair in atom_pair_list:
            if len(keep) == 0:
                keep.append(atom_pair)
            else:
                append = False
                remove = {}
                for compare in keep:
                    if atom_pair["acc index"] == compare["acc index"]:
                        if float(atom_pair["ang"]) > float(compare["ang"]):
                            remove = compare
                            append = True
                    else:
                        append = True
                if append:
                    keep.append(atom_pair)
                    if remove:
                        keep.remove(remove)
        for atom_pair in keep:
            writer.writerow([atom_pair["ang"], atom_pair["dist"]])

with open(f"acceptor_of_int_hbond.csv", "w") as csv_file:
    writer = csv.writer(csv_file)
    for acceptor in h_bond_to_acc_nonredundant.keys():
        for atom_pair_list in h_bond_to_acc_nonredundant[acceptor].values():
            for atom_pair in atom_pair_list:
                if atom_pair["vertex"] == "hydrogen":
                    writer.writerow([atom_pair["ang"], atom_pair["dist"]])

# count = 0
# values = []
for atom_pair_list in h_bond_from_don_nonredundant.values():
    keep = []
    for atom_pair in atom_pair_list:
        if len(keep) == 0:
            keep.append(atom_pair)
        else:
            append = False
            remove = {}
            for compare in keep:
                if atom_pair["acc index"] == compare["acc index"]:
                    if float(atom_pair["ang"]) > float(compare["ang"]):
                        remove = compare
                        append = True
                else:
                    append = True
            if append:
                keep.append(atom_pair)
                if remove:
                    keep.remove(remove)
    for atom_pair in keep:
        print(atom_pair)
    print()
    # values.append(len(x))
# print(len(nuc_data.keys()))
