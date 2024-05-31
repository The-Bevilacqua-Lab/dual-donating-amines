"""
This module will eventually do something.
"""

import pandas as pd
import csv
import const

# Set the H-bonding criteria.
H_DIST_MAX = 2.5
H_ANG_TOL = 40.0

# Construct a list of tuples containing atom and residue names that describe atoms capable of both donating and
# accepting an H-bond.
don_acc_atoms = []
for residue in const.RESIDUE_LIBRARY:
    for donor in residue['don']:
        for acceptor in residue['acc']:
            if donor[0] == acceptor[0]:
                don_acc_atoms.append((residue['res'], donor[0]))

# Extract the data from the h_bond csv file and remove redundant lines.
unique_col_comb = ["don_index", "acc_index", "h_name"]
# noinspection PyTypeChecker
h_bond_data = (pd.read_csv(snakemake.input.h_bond, comment="#", keep_default_na=False,
                           na_values={"h_acc_distance": "NaN", "h_angle": "NaN", "h_dihedral": "NaN"},
                           dtype={"don_chain": "object", "acc_chain": "object"})
               .drop_duplicates(subset=unique_col_comb))

# Extract the data from the nuc csv file.
nuc_data = pd.read_csv(snakemake.input.nuc, comment="#", na_filter=False, dtype={"chain": "object"})

# Create a dictionary based on the file containing information on nucleobase b-factors.
b_factor_grp = ["resn", "resi", "chain", "subset"]
b_factor_data = pd.read_csv(snakemake.input.b_factor, comment="#", na_filter=False, dtype={"chain": "object"})
b_factor_data["mean"] = b_factor_data.groupby(b_factor_grp)["b-factor"].transform("mean")
b_factor_data = b_factor_data.drop_duplicates(subset=b_factor_grp).drop(columns=["index", "name", "b-factor"])

# # Extract the data from the equivalence class file.
# eq_class_members_data = (pd.read_csv(snakemake.input.eq_class_members_data, header=None, comment="#", skiprows=[3],
#                                      na_filter=False, dtype="object")
#                          .set_index(pd.Series(["pdb_id", "model", "chain"])).transpose())
#
# # Identify atom pairs that meet the H-bond criteria and include a donor of interest.
# don_h_bonds = (h_bond_data[(h_bond_data["h_acc_distance"] <= H_DIST_MAX) & (h_bond_data["h_angle"] >= 180.0 - H_ANG_TOL)]
#               .merge(pd.DataFrame(const.DONORS_OF_INTEREST, columns=["don_resn", "don_name"]), how='inner')
#               .merge(eq_class_members_data, left_on="don_chain", right_on="chain", how='inner')
#               .drop(columns=["pdb_id", "model", "chain"]))
#
# # Identify atom pairs that meet the H-bond criteria and include a protonated donor of interest. Values in the _merge
# # column matching left_only indicate acceptors that cannot typically also donate an H-bond.
# prot_don_h_bonds_unfiltered = (h_bond_data[(h_bond_data["h_acc_distance"] <= H_DIST_MAX) & (h_bond_data["h_angle"] >= 180.0 - H_ANG_TOL)]
#                               .merge(pd.DataFrame(const.PROT_DONORS_OF_INTEREST, columns=["don_resn", "don_name"]),
#                                      how='inner')
#                               .merge(pd.DataFrame(don_acc_atoms, columns=["acc_resn", "acc_name"]), how='left',
#                                      indicator=True))
#
# # Filter prot_don_h_bonds_unfiltered such that only acceptors that cannot typically also donate an H-bond are included.
# # Also remove side chain acceptor atoms from ASN and GLN residues due to ambiguity in side chain conformation.
# prot_don_h_bonds = (prot_don_h_bonds_unfiltered[(prot_don_h_bonds_unfiltered["_merge"] == "left_only") &
#                                               ~((prot_don_h_bonds_unfiltered["acc_resn"] == "ASN") &
#                                                 (prot_don_h_bonds_unfiltered["acc_name"] == "OD1")) &
#                                               ~((prot_don_h_bonds_unfiltered["acc_resn"] == "GLN") &
#                                                 (prot_don_h_bonds_unfiltered["acc_name"] == "OE1"))]
#                    .drop(columns="_merge")
#                    .merge(eq_class_members_data, left_on="don_chain", right_on="chain", how='inner')
#                    .drop(columns=["pdb_id", "model", "chain"]))
#
# # Identify atom pairs that meet the H-bond criteria and include an acceptor of interest.
# acc_h_bonds = (h_bond_data[((h_bond_data["h_acc_distance"].notna()) & (h_bond_data["h_acc_distance"] <= H_DIST_MAX) &
#                           (h_bond_data["h_angle"] >= 180.0 - H_ANG_TOL)) |
#                          ((h_bond_data["h_acc_distance"].isna()) & (h_bond_data["don_acc_distance"] <= DON_DIST_MAX) &
#                                     (h_bond_data["don_angle"] >= DON_ANG_MIN) &
#                                     (h_bond_data["don_angle"] <= DON_ANG_MAX))]
#               .merge(pd.DataFrame(const.ACCEPTORS_OF_INTEREST, columns=["acc_resn", "acc_name"]), how='inner')
#               .merge(eq_class_members_data, left_on="acc_chain", right_on="chain", how='inner')
#               .drop(columns=["pdb_id", "model", "chain"]))
#
# # Create a new list of H-bonding atom pairs involving the donors of interest with redundancy removed such that if
# # multiple atom pairs involve the same donor and an acceptor belonging to different side chain conformations of the same
# # ASN, GLN, or HIS residue, only the atom pairs with the original conformation are kept. There could be situations where
# # one hydrogen of an exocyclic amine donates to the original side chain conformation of the acceptor residue while the
# # other hydrogen of the same exocyclic amine donates to the rotated side chain conformation of the acceptor residue. For
# # simplicity, H-bonds involving the original side chain conformations of acceptor residues are prioritized over
# # situations where a dual H-bonding exocyclic amine could be inferred by discarding the atom pair involving the original
# # side chain conformation and keeping the atom pair involving the rotated side chain conformation. Additionally, create
# # a new list of H-bonding atom pairs involving the protonated donors of interest with redundancy removed.
# acc_grp = ["don_index", "acc_resn", "acc_resi", "acc_chain"]
# don_h_bonds_nr = (don_h_bonds[don_h_bonds.groupby(acc_grp)["rotated_side_chain"]
#                  .transform(lambda grp: grp.str.fullmatch("NaN") if any(grp.str.fullmatch("NaN")) else True)])
# prot_don_h_bonds_nr = (prot_don_h_bonds[prot_don_h_bonds.groupby(acc_grp)["rotated_side_chain"]
#                       .transform(lambda grp: grp.str.fullmatch("NaN") if any(grp.str.fullmatch("NaN")) else True)])
#
# # Create a new list of H-bonding atom pairs involving the acceptors of interest with redundancy removed such that if
# # multiple atom pairs involve the same acceptor and a donor belonging to different side chain conformations of the same
# # ASN, GLN, or HIS residue, only the atom pairs with the original conformation are kept. Redundancy pertaining to
# # donation by TYR OH atoms is also removed. The TYR OH hydrogen can occupy two different conformations. If both
# # conformations are capable of forming an H-bond with the acceptor of interest, only the conformation with the shortest
# # hydrogen-acceptor distance is kept. Additionally, create a new list of H-bonding atom pairs involving the deprotonated
# # acceptors of interest with redundancy removed.
# don_grp = ["don_resn", "don_resi", "don_chain", "acc_index"]
# acc_h_bonds_nr = (acc_h_bonds[(((acc_h_bonds["don_resn"] != "TYR") &
#                               (acc_h_bonds.groupby(don_grp)["rotated_side_chain"]
#                                .transform(lambda grp: grp.str.fullmatch("NaN")
#                                           if any(grp.str.fullmatch("NaN")) else True))) |
#                              ((acc_h_bonds["don_resn"] == "TYR") &
#                               (acc_h_bonds.groupby(don_grp)["h_acc_distance"]
#                                .transform(lambda grp: [mem == grp.min() for mem in grp]))))])
#
# # Write dataframes to csv files.
# nuc_data.assign(eq_class_members=snakemake.wildcards.eq_class_members).to_csv(snakemake.output.nuc_data, index=False)
# don_h_bonds_nr.assign(eq_class_members=snakemake.wildcards.eq_class_members).to_csv(snakemake.output.don_h_bonds_nr, index=False)
# prot_don_h_bonds_nr.assign(eq_class_members=snakemake.wildcards.eq_class_members).to_csv(snakemake.output.prot_don_h_bonds_nr, index=False)
# acc_h_bonds_nr.assign(eq_class_members=snakemake.wildcards.eq_class_members).to_csv(snakemake.output.acc_h_bonds_nr, index=False)
# b_factor_data.assign(eq_class_members=snakemake.wildcards.eq_class_members).to_csv(snakemake.output.b_factor_data, index=False)
# h_bond_data.assign(eq_class_members=snakemake.wildcards.eq_class_members).to_csv(snakemake.output.h_bond_data, index=False)
