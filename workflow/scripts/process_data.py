"""
This module will eventually do something.
"""

import pandas as pd
import csv
import const

# Set the H-bonding criteria.
H_DIST_MAX = 2.5
H_ANG_TOL = 40.0
DON_DIST_MAX = 3.3
DON_ANG_MIN = 100.0
DON_ANG_MAX = 145.0

# Construct a list of tuples containing atom and residue names that describe atoms capable of both donating and
# accepting an H-bond.
don_acc_atoms = []
for residue in const.RESIDUE_LIBRARY:
    for donor in residue['don']:
        for acceptor in residue['acc']:
            if donor[0] == acceptor[0]:
                don_acc_atoms.append((residue['res'], donor[0]))

# Extract the data from the hbond csv file and remove redundant lines.
hbond_col_names = ["don_index", "don_name", "don_resn", "don_resi", "don_chain", "acc_index", "acc_name",
                   "acc_resn", "acc_resi", "acc_chain", "don_acc_distance", "don_angle", "acc_angle", "h_acc_distance",
                   "h_angle", "h_dihedral", "h_name", "rotated_side_chain"]
unique_col_comb = ["don_index", "acc_index", "h_name", "rotated_side_chain"]
# noinspection PyTypeChecker
hbond_data = (pd.read_csv(snakemake.input.hbond_data, names=hbond_col_names, comment="#", keep_default_na=False,
                          na_values={"h_acc_distance": "NaN", "h_angle": "NaN", "h_dihedral": "NaN"},
                          dtype={"don_chain": "object", "acc_chain": "object"})
              .drop_duplicates(subset=unique_col_comb))

# Extract the data from the nuc csv file.
nuc_col_names = ["index", "atom_name", "resn", "resi", "chain"]
nuc_data = pd.read_csv(snakemake.input.nuc_data, names=nuc_col_names, comment="#", na_filter=False,
                       dtype={"chain": "object"})

# Create a dictionary based on the file containing information on nucleobase b-factors.
with open(snakemake.input.b_factor_data, mode='r') as read_file:
    b_factors_dict = {
        "resn": [],
        "resi": [],
        "chain": [],
        "b_atom": []
    }
    for line in csv.reader(read_file):
        if line[0][0] != "#":
            for i in range(len(line) - 3):
                b_factors_dict["resn"].append(line[0])
                b_factors_dict["resi"].append(line[1])
                b_factors_dict["chain"].append(line[2])
                b_factors_dict["b_atom"].append(line[i + 3])
b_factor_data = pd.DataFrame(b_factors_dict).astype({"b_atom": "float64"})
b_factor_data["b_mean"] = b_factor_data.groupby(["resn", "resi", "chain"])["b_atom"].transform("mean")
b_factor_data = b_factor_data.drop_duplicates(subset=["resn", "resi", "chain"]).drop(columns=["b_atom"])

# Extract the data from the equivalence class file.
eq_class_members_data = (pd.read_csv(snakemake.input.eq_class_members_data, header=None, comment="#", skiprows=[3],
                                     na_filter=False, dtype="object")
                         .set_index(pd.Series(["pdb_id", "model", "chain"])).transpose())

# Identify atom pairs that meet the H-bond criteria and include a donor of interest.
don_hbonds = (hbond_data[(hbond_data["h_acc_distance"] <= H_DIST_MAX) & (hbond_data["h_angle"] >= 180.0 - H_ANG_TOL)]
              .merge(pd.DataFrame(const.DONORS_OF_INTEREST, columns=["don_resn", "don_name"]), how='inner')
              .merge(eq_class_members_data, left_on="don_chain", right_on="chain", how='inner')
              .drop(columns=["pdb_id", "model", "chain"]))

# Identify atom pairs that meet the H-bond criteria and include a protonated donor of interest. Values in the _merge
# column matching left_only indicate acceptors that cannot typically also donate an H-bond.
prot_don_hbonds_unfiltered = (hbond_data[(hbond_data["h_acc_distance"] <= H_DIST_MAX) & (hbond_data["h_angle"] >= 180.0 - H_ANG_TOL)]
                              .merge(pd.DataFrame(const.PROT_DONORS_OF_INTEREST, columns=["don_resn", "don_name"]),
                                     how='inner')
                              .merge(pd.DataFrame(don_acc_atoms, columns=["acc_resn", "acc_name"]), how='left',
                                     indicator=True))

# Filter prot_don_hbonds_unfiltered such that only acceptors that cannot typically also donate an H-bond are included.
# Also remove side chain acceptor atoms from ASN and GLN residues due to ambiguity in side chain conformation.
prot_don_hbonds = (prot_don_hbonds_unfiltered[(prot_don_hbonds_unfiltered["_merge"] == "left_only") &
                                              ~((prot_don_hbonds_unfiltered["acc_resn"] == "ASN") &
                                                (prot_don_hbonds_unfiltered["acc_name"] == "OD1")) &
                                              ~((prot_don_hbonds_unfiltered["acc_resn"] == "GLN") &
                                                (prot_don_hbonds_unfiltered["acc_name"] == "OE1"))]
                   .drop(columns="_merge")
                   .merge(eq_class_members_data, left_on="don_chain", right_on="chain", how='inner')
                   .drop(columns=["pdb_id", "model", "chain"]))

# Identify atom pairs that meet the H-bond criteria and include an acceptor of interest.
acc_hbonds = (hbond_data[((hbond_data["h_acc_distance"].notna()) & (hbond_data["h_acc_distance"] <= H_DIST_MAX) &
                          (hbond_data["h_angle"] >= 180.0 - H_ANG_TOL)) |
                         ((hbond_data["h_acc_distance"].isna()) & (hbond_data["don_acc_distance"] <= DON_DIST_MAX) &
                                    (hbond_data["don_angle"] >= DON_ANG_MIN) &
                                    (hbond_data["don_angle"] <= DON_ANG_MAX))]
              .merge(pd.DataFrame(const.ACCEPTORS_OF_INTEREST, columns=["acc_resn", "acc_name"]), how='inner')
              .merge(eq_class_members_data, left_on="acc_chain", right_on="chain", how='inner')
              .drop(columns=["pdb_id", "model", "chain"]))

# Create a new list of H-bonding atom pairs involving the donors of interest with redundancy removed such that if
# multiple atom pairs involve the same donor and an acceptor belonging to different side chain conformations of the same
# ASN, GLN, or HIS residue, only the atom pairs with the original conformation are kept. There could be situations where
# one hydrogen of an exocyclic amine donates to the original side chain conformation of the acceptor residue while the
# other hydrogen of the same exocyclic amine donates to the rotated side chain conformation of the acceptor residue. For
# simplicity, H-bonds involving the original side chain conformations of acceptor residues are prioritized over
# situations where a dual H-bonding exocyclic amine could be inferred by discarding the atom pair involving the original
# side chain conformation and keeping the atom pair involving the rotated side chain conformation. Additionally, create
# a new list of H-bonding atom pairs involving the protonated donors of interest with redundancy removed.
acc_grp = ["don_index", "acc_resn", "acc_resi", "acc_chain"]
don_hbonds_nr = (don_hbonds[don_hbonds.groupby(acc_grp)["rotated_side_chain"]
                 .transform(lambda grp: grp.str.fullmatch("NaN") if any(grp.str.fullmatch("NaN")) else True)])
prot_don_hbonds_nr = (prot_don_hbonds[prot_don_hbonds.groupby(acc_grp)["rotated_side_chain"]
                      .transform(lambda grp: grp.str.fullmatch("NaN") if any(grp.str.fullmatch("NaN")) else True)])

# Create a new list of H-bonding atom pairs involving the acceptors of interest with redundancy removed such that if
# multiple atom pairs involve the same acceptor and a donor belonging to different side chain conformations of the same
# ASN, GLN, or HIS residue, only the atom pairs with the original conformation are kept. Redundancy pertaining to
# donation by TYR OH atoms is also removed. The TYR OH hydrogen can occupy two different conformations. If both
# conformations are capable of forming an H-bond with the acceptor of interest, only the conformation with the shortest
# hydrogen-acceptor distance is kept. Additionally, create a new list of H-bonding atom pairs involving the deprotonated
# acceptors of interest with redundancy removed.
don_grp = ["don_resn", "don_resi", "don_chain", "acc_index"]
acc_hbonds_nr = (acc_hbonds[(((acc_hbonds["don_resn"] != "TYR") &
                              (acc_hbonds.groupby(don_grp)["rotated_side_chain"]
                               .transform(lambda grp: grp.str.fullmatch("NaN")
                                          if any(grp.str.fullmatch("NaN")) else True))) |
                             ((acc_hbonds["don_resn"] == "TYR") &
                              (acc_hbonds.groupby(don_grp)["h_acc_distance"]
                               .transform(lambda grp: [mem == grp.min() for mem in grp]))))])

# Write dataframes to csv files.
nuc_data.assign(eq_class_members=snakemake.wildcards.eq_class_members).to_csv(snakemake.output.nuc_data, index=False)
don_hbonds_nr.assign(eq_class_members=snakemake.wildcards.eq_class_members).to_csv(snakemake.output.don_hbonds_nr, index=False)
prot_don_hbonds_nr.assign(eq_class_members=snakemake.wildcards.eq_class_members).to_csv(snakemake.output.prot_don_hbonds_nr, index=False)
acc_hbonds_nr.assign(eq_class_members=snakemake.wildcards.eq_class_members).to_csv(snakemake.output.acc_hbonds_nr, index=False)
b_factor_data.assign(eq_class_members=snakemake.wildcards.eq_class_members).to_csv(snakemake.output.b_factor_data, index=False)
hbond_data.assign(eq_class_members=snakemake.wildcards.eq_class_members).to_csv(snakemake.output.hbond_data, index=False)

# TODO prot donor should not donate to the same atom the exo amine is donating to (e.g., C788 A5 in 6xu8)
# TODO keep in mind ["H04", "A", "4212", "L5", "OD1", "ASN", "3", "LT"] in 8GLP
# TODO why is NR_3.0_65061.22_hbond_test_data/LA/LA/ARG`3/NH2 missing a proton?
