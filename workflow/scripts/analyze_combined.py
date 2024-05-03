"""
This script will eventually do something.
"""

import os
import pandas as pd
import const

# the working directory should house the combined/ folder
working_dir = os.getcwd()

# prepare a list of acceptor residue names and atom names that have substantially greater negative charge
neg_acc_resn = []
neg_acc_name = []
for residue in const.RESIDUE_LIBRARY:
    for acceptor in residue['acc']:
        if acceptor[3] < 0:
            neg_acc_resn.append(residue['res'])
            neg_acc_name.append(acceptor[0])
        elif residue['res'] == 'ASP' and acceptor[0] == 'OD1':
            neg_acc_resn.append(residue['res'])
            neg_acc_name.append(acceptor[0])
        elif residue['res'] == 'GLU' and acceptor[0] == 'OE1':
            neg_acc_resn.append(residue['res'])
            neg_acc_name.append(acceptor[0])
        elif residue['res'] in ['A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT'] and acceptor[0] == 'OP1':
            neg_acc_resn.append(residue['res'])
            neg_acc_name.append(acceptor[0])

# prepare a dataframe of acceptors that have substantially negative charges
neg_acceptors = pd.DataFrame({
    "acc_resn": neg_acc_resn,
    "acc_name": neg_acc_name
})

# create the combined dataframes
don_hbonds_nr = pd.read_csv(working_dir + "/combined/don_hbonds_nr_c.csv", na_filter=False,
                            dtype={"don_resi": "object", "acc_resi": "object"})
prot_don_hbonds_nr = pd.read_csv(working_dir + "/combined/prot_don_hbonds_nr_c.csv", na_filter=False,
                                 dtype={"don_resi": "object", "acc_resi": "object"})
acc_hbonds_nr = pd.read_csv(working_dir + "/combined/acc_hbonds_nr_c.csv", na_filter=False,
                            dtype={"don_resi": "object", "acc_resi": "object"})
deprot_acc_hbonds_nr = pd.read_csv(working_dir + "/combined/deprot_acc_hbonds_nr_c.csv", na_filter=False,
                                   dtype={"don_resi": "object", "acc_resi": "object"})
nuc_data_raw = pd.read_csv(working_dir + "/combined/nuc_data_c.csv", na_filter=False, dtype={"resi": "object"})
b_factor_data = pd.read_csv(working_dir + "/combined/b_factor_data_c.csv", na_filter=False, dtype={"resi": "object"})
hbond_data = pd.read_csv(working_dir + "/combined/hbond_data_c.csv", na_filter=False,
                         dtype={"don_resi": "object", "acc_resi": "object"})

# Create a new dataframe of nucleobases where the mean of the b-factors are below 79.
nuc_data = nuc_data_raw.merge(b_factor_data[b_factor_data["b_mean"] < 79.0], on=["resn", "resi", "chain",
                                                                                 "eq_class_members"],
                              how='inner').drop(columns=["b_mean"])

# Prepare data on donor-acceptor pairs that can be used to look at the distribution of h-bonding geometry measurements.
acc_grp = ["don_index", "acc_resn", "acc_resi", "acc_chain", "eq_class_members"]
don_grp = ["acc_index", "don_resn", "don_resi", "don_chain", "eq_class_members"]
# Remove redundancy in acceptor residues for each donor.
hbonds_geom_nr_acc_only = (hbond_data[hbond_data.groupby(acc_grp)["rotated_side_chain"]
                           .transform(
                               lambda grp: grp.str.fullmatch("none") if any(grp.str.fullmatch("none")) else True)])
# Remove redundancy in donor residues for each acceptor.
hbonds_geom_nr = (hbonds_geom_nr_acc_only[hbonds_geom_nr_acc_only.groupby(don_grp)["rotated_side_chain"]
                  .transform(lambda grp: grp.str.fullmatch("none") if any(grp.str.fullmatch("none")) else True)])

# Write data on donor-acceptor pairs to csv files.
don_acc_grp = ["don_index", "acc_index", "eq_class_members"]
(hbonds_geom_nr[
     # only include hydrogens with smaller D-H...A distances
     (hbonds_geom_nr.groupby(don_acc_grp)["h_acc_distance"]
      .transform(lambda grp: [mem == grp.min() for mem in grp])) &
     # only consider non-rotatable donors
     (hbonds_geom_nr["h_acc_distance"] != "NaN")]
 .to_csv("plots/geom/hbonds_geom_nr_non_rot.csv", index=False, columns=["h_acc_distance", "h_angle"]))
(hbonds_geom_nr[
     # only include hydrogens with smaller D-H...A distances
     (hbonds_geom_nr.groupby(don_acc_grp)["h_acc_distance"]
      .transform(lambda grp: [mem == grp.min() for mem in grp])) &
     # only consider non-rotatable donors
     (hbonds_geom_nr["h_acc_distance"] != "NaN") &
     # do not consider A(N6)-U(O4), C(N4)-G(O6), G(N2)-C(O2), G(N2)-C(N3), U(N3)-A(N1), G(N1)-C(N3), and G(N1)-C(O2)
     # atom pairs
     (~hbonds_geom_nr[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["A", "N6", "U", "O4"])
      .all(axis='columns')) &
     (~hbonds_geom_nr[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["C", "N4", "G", "O6"])
      .all(axis='columns')) &
     (~hbonds_geom_nr[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["G", "N2", "C", "O2"])
      .all(axis='columns')) &
     (~hbonds_geom_nr[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["G", "N2", "C", "N3"])
      .all(axis='columns')) &
     (~hbonds_geom_nr[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["U", "N3", "A", "N1"])
      .all(axis='columns')) &
     (~hbonds_geom_nr[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["G", "N1", "C", "N3"])
      .all(axis='columns')) &
     (~hbonds_geom_nr[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["G", "N1", "C", "O2"])
      .all(axis='columns'))]
 .to_csv("plots/geom/hbonds_geom_nr_non_rot_filtered.csv", index=False, columns=["h_acc_distance", "h_angle"]))
don_acc_grp = ["don_index", "acc_index", "eq_class_members"]
(hbonds_geom_nr[
     # only include pairs with smaller D-A distances
     (hbonds_geom_nr.groupby(don_acc_grp)["don_acc_distance"]
      .transform(lambda grp: [mem == grp.min() for mem in grp])) &
     # only consider non-rotatable donors
     (hbonds_geom_nr["h_acc_distance"] == "NaN")]
 .to_csv("plots/geom/hbonds_geom_nr_rot.csv", index=False, columns=["don_acc_distance", "don_angle"]))

# Identify nucleobases that donate either one or at least two H-bonds via their exocyclic amines. The H-bonds from the
# latter nucleobases must involve both exocyclic amine hydrogens and at least two different acceptors.
single_don_hbonds = don_hbonds_nr[
    (don_hbonds_nr.groupby(["don_index", "eq_class_members"])["h_name"].transform("nunique") == 1) |
    (don_hbonds_nr.groupby(["don_index", "eq_class_members"])["acc_index"].transform("nunique") == 1)]
dual_don_hbonds = don_hbonds_nr[
    (don_hbonds_nr.groupby(["don_index", "eq_class_members"])["h_name"].transform("nunique") == 2) &
    (don_hbonds_nr.groupby(["don_index", "eq_class_members"])["acc_index"].transform("nunique") >= 2)]

# Prepare dataframes of exocyclic amines for each A, C, and G residue.
a_exo_amines = nuc_data[nuc_data["atom_name"] == "N6"]

# Prepare dataframes of single and dual H-bonding A(N6) with H-bond information.
a_n6_single_hbond = (single_don_hbonds[single_don_hbonds["don_resn"] == "A"]
                     .merge(a_exo_amines,
                            left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                            right_on=["resn", "resi", "chain", "eq_class_members"],
                            how='inner').drop(columns=["index", "atom_name", "resn", "resi", "chain"]))
a_n6_dual_hbond = (dual_don_hbonds[dual_don_hbonds["don_resn"] == "A"]
                   .merge(a_exo_amines,
                          left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                          right_on=["resn", "resi", "chain", "eq_class_members"],
                          how='inner').drop(columns=["index", "atom_name", "resn", "resi", "chain"]))

# Prepare dataframes of no, single, and dual H-bonding A(N6) with just residue information.
a_n6_single_res = (a_n6_single_hbond.loc[:, ["eq_class_members", "don_resn", "don_resi", "don_chain"]]
                   .drop_duplicates()).rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"})
a_n6_dual_res = (a_n6_dual_hbond.loc[:, ["eq_class_members", "don_resn", "don_resi", "don_chain"]]
                 .drop_duplicates()).rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"})
a_n6_no_res = (pd.concat([a_n6_single_res,
                         a_n6_dual_res,
                         a_exo_amines.loc[:, ["eq_class_members", "resn", "resi", "chain"]]])
               .drop_duplicates(keep=False))

# Prepare a list of single and dual H-bonding A residues that are involved in canonical AU base pairs.
au_bp_single_res = (a_n6_single_hbond[a_n6_single_hbond[["acc_name", "acc_resn"]].eq(["O4", "U"]).all(axis='columns')]
                    .merge(acc_hbonds_nr[acc_hbonds_nr[["don_name", "don_resn", "acc_name"]].eq(["N3", "U", "N1"])
                           .all(axis='columns')],
                           left_on=["don_resn", "don_resi", "don_chain", "acc_resn", "acc_resi", "acc_chain",
                                    "eq_class_members"],
                           right_on=["acc_resn", "acc_resi", "acc_chain", "don_resn", "don_resi", "don_chain",
                                     "eq_class_members"],
                           how='inner').loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                    .drop_duplicates()
                    .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"}))
au_bp_dual_res = (a_n6_dual_hbond[a_n6_dual_hbond[["acc_name", "acc_resn"]].eq(["O4", "U"]).all(axis='columns')]
                  .merge(acc_hbonds_nr[acc_hbonds_nr[["don_name", "don_resn", "acc_name"]].eq(["N3", "U", "N1"])
                         .all(axis='columns')],
                         left_on=["don_resn", "don_resi", "don_chain", "acc_resn", "acc_resi", "acc_chain",
                                  "eq_class_members"],
                         right_on=["acc_resn", "acc_resi", "acc_chain", "don_resn", "don_resi", "don_chain",
                                   "eq_class_members"],
                         how='inner').loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                  .drop_duplicates()
                  .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"}))

# Write data regarding single and dual H-bonding A residues that are involved in canonical AU base pairs.
au_bp = pd.DataFrame({
    "Type": ["Single", "Dual"],
    "Occurrence": [au_bp_single_res.shape[0],
                   au_bp_dual_res.shape[0]],
    "Total": [a_n6_single_res.shape[0],
              a_n6_dual_res.shape[0]]
})
au_bp.to_csv("plots/frequency/au_bp.csv", index=False)

# Prepare dataframes with just residue information of A residues that accept H-bonds at their N1, N3, or N7.
a_n1_acc_res = (acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]].eq(["N1", "A"]).all(axis='columns')]
                .loc[:, ["eq_class_members", "acc_resn", "acc_resi", "acc_chain"]].drop_duplicates()
                .rename(columns={"acc_resn": "resn", "acc_resi": "resi", "acc_chain": "chain"}))
a_n3_acc_res = (acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]].eq(["N3", "A"]).all(axis='columns')]
                .loc[:, ["eq_class_members", "acc_resn", "acc_resi", "acc_chain"]].drop_duplicates()
                .rename(columns={"acc_resn": "resn", "acc_resi": "resi", "acc_chain": "chain"}))
a_n7_acc_res = (acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]].eq(["N7", "A"]).all(axis='columns')]
                .loc[:, ["eq_class_members", "acc_resn", "acc_resi", "acc_chain"]].drop_duplicates()
                .rename(columns={"acc_resn": "resn", "acc_resi": "resi", "acc_chain": "chain"}))

# Prepare dataframes of no, single, and dual H-bonding A(N6) that also accept at the N1 with just residue information.
a_n6_no_a_n1_acc_res = a_n6_no_res.merge(a_n1_acc_res, how='inner')
a_n6_single_a_n1_acc_res = a_n6_single_res.merge(a_n1_acc_res, how='inner')
a_n6_dual_a_n1_acc_res = a_n6_dual_res.merge(a_n1_acc_res, how='inner')

# Prepare dataframes of no, single, and dual H-bonding A(N6) that also accept at the N3 with just residue information.
a_n6_no_a_n3_acc_res = a_n6_no_res.merge(a_n3_acc_res, how='inner')
a_n6_single_a_n3_acc_res = a_n6_single_res.merge(a_n3_acc_res, how='inner')
a_n6_dual_a_n3_acc_res = a_n6_dual_res.merge(a_n3_acc_res, how='inner')

# Prepare dataframes of no, single, and dual H-bonding A(N6) that also accept at the N7 with just residue information.
a_n6_no_a_n7_acc_res = a_n6_no_res.merge(a_n7_acc_res, how='inner')
a_n6_single_a_n7_acc_res = a_n6_single_res.merge(a_n7_acc_res, how='inner')
a_n6_dual_a_n7_acc_res = a_n6_dual_res.merge(a_n7_acc_res, how='inner')

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at N1 and which include H-bond information
# related to the N6 H-bond donation.
a_n6_single_a_n1_acc_hbond_from_n6 = (a_n6_single_a_n1_acc_res
                                      .rename(columns={"resn": "don_resn", "resi": "don_resi", "chain": "don_chain"})
                                      .merge(a_n6_single_hbond, how='inner'))
a_n6_dual_a_n1_acc_hbond_from_n6 = (a_n6_dual_a_n1_acc_res
                                    .rename(columns={"resn": "don_resn", "resi": "don_resi", "chain": "don_chain"})
                                    .merge(a_n6_dual_hbond, how='inner'))

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at N3 and which include H-bond information
# related to the N6 H-bond donation.
a_n6_single_a_n3_acc_hbond_from_n6 = (a_n6_single_a_n3_acc_res
                                      .rename(columns={"resn": "don_resn", "resi": "don_resi", "chain": "don_chain"})
                                      .merge(a_n6_single_hbond, how='inner'))
a_n6_dual_a_n3_acc_hbond_from_n6 = (a_n6_dual_a_n3_acc_res
                                    .rename(columns={"resn": "don_resn", "resi": "don_resi", "chain": "don_chain"})
                                    .merge(a_n6_dual_hbond, how='inner'))

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at N7 and which include H-bond information
# related to the N6 H-bond donation.
a_n6_single_a_n7_acc_hbond_from_n6 = (a_n6_single_a_n7_acc_res
                                      .rename(columns={"resn": "don_resn", "resi": "don_resi", "chain": "don_chain"})
                                      .merge(a_n6_single_hbond, how='inner'))
a_n6_dual_a_n7_acc_hbond_from_n6 = (a_n6_dual_a_n7_acc_res
                                    .rename(columns={"resn": "don_resn", "resi": "don_resi", "chain": "don_chain"})
                                    .merge(a_n6_dual_hbond, how='inner'))

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at N1 and which include H-bond information
# related to the N1 H-bond acceptation.
a_n6_single_a_n1_acc_hbond_to_n1 = (a_n6_single_a_n1_acc_res
                                    .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                    .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]].eq(["N1", "A"])
                                           .all(axis='columns')], how='inner'))
a_n6_dual_a_n1_acc_hbond_to_n1 = (a_n6_dual_a_n1_acc_res
                                  .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                  .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]].eq(["N1", "A"])
                                         .all(axis='columns')], how='inner'))

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at N3 and which include H-bond information
# related to the N3 H-bond acceptation.
a_n6_single_a_n3_acc_hbond_to_n3 = (a_n6_single_a_n3_acc_res
                                    .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                    .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]].eq(["N3", "A"])
                                           .all(axis='columns')], how='inner'))
a_n6_dual_a_n3_acc_hbond_to_n3 = (a_n6_dual_a_n3_acc_res
                                  .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                  .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]].eq(["N3", "A"])
                                         .all(axis='columns')], how='inner'))

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at N7 and which include H-bond information
# related to the N7 H-bond acceptation.
a_n6_single_a_n7_acc_hbond_to_n7 = (a_n6_single_a_n7_acc_res
                                    .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                    .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]].eq(["N7", "A"])
                                           .all(axis='columns')], how='inner'))
a_n6_dual_a_n7_acc_hbond_to_n7 = (a_n6_dual_a_n7_acc_res
                                  .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                  .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]].eq(["N7", "A"])
                                         .all(axis='columns')], how='inner'))

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at N1 and which include H-bond information
# related to the N6 H-bond donation and the N1 H-bond acceptation.
a_n6_single_a_n1_acc_match = (a_n6_single_a_n1_acc_hbond_from_n6
                              .merge(a_n6_single_a_n1_acc_hbond_to_n1,
                                     left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                     right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"], how='inner'))
a_n6_dual_a_n1_acc_match = (a_n6_dual_a_n1_acc_hbond_from_n6
                            .merge(a_n6_dual_a_n1_acc_hbond_to_n1,
                                   left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                   right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"], how='inner'))

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at N3 and which include H-bond information
# related to the N6 H-bond donation and the N3 H-bond acceptation.
a_n6_single_a_n3_acc_match = (a_n6_single_a_n3_acc_hbond_from_n6
                              .merge(a_n6_single_a_n3_acc_hbond_to_n3,
                                     left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                     right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"], how='inner'))
a_n6_dual_a_n3_acc_match = (a_n6_dual_a_n3_acc_hbond_from_n6
                            .merge(a_n6_dual_a_n3_acc_hbond_to_n3,
                                   left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                   right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"], how='inner'))

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at N7 and which include H-bond information
# related to the N6 H-bond donation and the N7 H-bond acceptation.
a_n6_single_a_n7_acc_match = (a_n6_single_a_n7_acc_hbond_from_n6
                              .merge(a_n6_single_a_n7_acc_hbond_to_n7,
                                     left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                     right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"], how='inner'))
a_n6_dual_a_n7_acc_match = (a_n6_dual_a_n7_acc_hbond_from_n6
                            .merge(a_n6_dual_a_n7_acc_hbond_to_n7,
                                   left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                   right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"], how='inner'))

# Prepare dataframes with just residue information of single and dual H-bonding A(N6) that also accepts an H-bond at the
# N1, N3, or N7. Additionally, there is overlap in the identity of the partner nucleobase, amino acid residue backbone,
# or amino acid residue side chain, where the N6 donates an H-bond and the N1, N3, or N7 accepts an H-bond from the same
# partner entity.
a_n6_single_a_n1_acc_ol = (a_n6_single_a_n1_acc_match[
    (((a_n6_single_a_n1_acc_match["acc_resn_x"] == a_n6_single_a_n1_acc_match["don_resn_y"]) &
      (a_n6_single_a_n1_acc_match["acc_resi_x"] == a_n6_single_a_n1_acc_match["don_resi_y"]) &
      (a_n6_single_a_n1_acc_match["acc_chain_x"] == a_n6_single_a_n1_acc_match["don_chain_y"]) &
      ~(a_n6_single_a_n1_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_single_a_n1_acc_match["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_single_a_n1_acc_match["acc_resn_x"] == a_n6_single_a_n1_acc_match["don_resn_y"]) &
      (a_n6_single_a_n1_acc_match["acc_resi_x"] == a_n6_single_a_n1_acc_match["don_resi_y"]) &
      (a_n6_single_a_n1_acc_match["acc_chain_x"] == a_n6_single_a_n1_acc_match["don_chain_y"]) &
      ~(a_n6_single_a_n1_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_single_a_n1_acc_match["don_name_y"].isin(["O2'"])) &
      (a_n6_single_a_n1_acc_match["acc_name_x"].isin(["O"])) &
      (a_n6_single_a_n1_acc_match["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())
a_n6_dual_a_n1_acc_ol = (a_n6_dual_a_n1_acc_match[
    (((a_n6_dual_a_n1_acc_match["acc_resn_x"] == a_n6_dual_a_n1_acc_match["don_resn_y"]) &
      (a_n6_dual_a_n1_acc_match["acc_resi_x"] == a_n6_dual_a_n1_acc_match["don_resi_y"]) &
      (a_n6_dual_a_n1_acc_match["acc_chain_x"] == a_n6_dual_a_n1_acc_match["don_chain_y"]) &
      ~(a_n6_dual_a_n1_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_dual_a_n1_acc_match["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_dual_a_n1_acc_match["acc_resn_x"] == a_n6_dual_a_n1_acc_match["don_resn_y"]) &
      (a_n6_dual_a_n1_acc_match["acc_resi_x"] == a_n6_dual_a_n1_acc_match["don_resi_y"]) &
      (a_n6_dual_a_n1_acc_match["acc_chain_x"] == a_n6_dual_a_n1_acc_match["don_chain_y"]) &
      ~(a_n6_dual_a_n1_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_dual_a_n1_acc_match["don_name_y"].isin(["O2'"])) &
      (a_n6_dual_a_n1_acc_match["acc_name_x"].isin(["O"])) &
      (a_n6_dual_a_n1_acc_match["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())
a_n6_single_a_n3_acc_ol = (a_n6_single_a_n3_acc_match[
    (((a_n6_single_a_n3_acc_match["acc_resn_x"] == a_n6_single_a_n3_acc_match["don_resn_y"]) &
      (a_n6_single_a_n3_acc_match["acc_resi_x"] == a_n6_single_a_n3_acc_match["don_resi_y"]) &
      (a_n6_single_a_n3_acc_match["acc_chain_x"] == a_n6_single_a_n3_acc_match["don_chain_y"]) &
      ~(a_n6_single_a_n3_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_single_a_n3_acc_match["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_single_a_n3_acc_match["acc_resn_x"] == a_n6_single_a_n3_acc_match["don_resn_y"]) &
      (a_n6_single_a_n3_acc_match["acc_resi_x"] == a_n6_single_a_n3_acc_match["don_resi_y"]) &
      (a_n6_single_a_n3_acc_match["acc_chain_x"] == a_n6_single_a_n3_acc_match["don_chain_y"]) &
      ~(a_n6_single_a_n3_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_single_a_n3_acc_match["don_name_y"].isin(["O2'"])) &
      (a_n6_single_a_n3_acc_match["acc_name_x"].isin(["O"])) &
      (a_n6_single_a_n3_acc_match["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())
a_n6_dual_a_n3_acc_ol = (a_n6_dual_a_n3_acc_match[
    (((a_n6_dual_a_n3_acc_match["acc_resn_x"] == a_n6_dual_a_n3_acc_match["don_resn_y"]) &
      (a_n6_dual_a_n3_acc_match["acc_resi_x"] == a_n6_dual_a_n3_acc_match["don_resi_y"]) &
      (a_n6_dual_a_n3_acc_match["acc_chain_x"] == a_n6_dual_a_n3_acc_match["don_chain_y"]) &
      ~(a_n6_dual_a_n3_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_dual_a_n3_acc_match["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_dual_a_n3_acc_match["acc_resn_x"] == a_n6_dual_a_n3_acc_match["don_resn_y"]) &
      (a_n6_dual_a_n3_acc_match["acc_resi_x"] == a_n6_dual_a_n3_acc_match["don_resi_y"]) &
      (a_n6_dual_a_n3_acc_match["acc_chain_x"] == a_n6_dual_a_n3_acc_match["don_chain_y"]) &
      ~(a_n6_dual_a_n3_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_dual_a_n3_acc_match["don_name_y"].isin(["O2'"])) &
      (a_n6_dual_a_n3_acc_match["acc_name_x"].isin(["O"])) &
      (a_n6_dual_a_n3_acc_match["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())
a_n6_single_a_n7_acc_ol = (a_n6_single_a_n7_acc_match[
    (((a_n6_single_a_n7_acc_match["acc_resn_x"] == a_n6_single_a_n7_acc_match["don_resn_y"]) &
      (a_n6_single_a_n7_acc_match["acc_resi_x"] == a_n6_single_a_n7_acc_match["don_resi_y"]) &
      (a_n6_single_a_n7_acc_match["acc_chain_x"] == a_n6_single_a_n7_acc_match["don_chain_y"]) &
      ~(a_n6_single_a_n7_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_single_a_n7_acc_match["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_single_a_n7_acc_match["acc_resn_x"] == a_n6_single_a_n7_acc_match["don_resn_y"]) &
      (a_n6_single_a_n7_acc_match["acc_resi_x"] == a_n6_single_a_n7_acc_match["don_resi_y"]) &
      (a_n6_single_a_n7_acc_match["acc_chain_x"] == a_n6_single_a_n7_acc_match["don_chain_y"]) &
      ~(a_n6_single_a_n7_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_single_a_n7_acc_match["don_name_y"].isin(["O2'"])) &
      (a_n6_single_a_n7_acc_match["acc_name_x"].isin(["O"])) &
      (a_n6_single_a_n7_acc_match["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())
a_n6_dual_a_n7_acc_ol = (a_n6_dual_a_n7_acc_match[
    (((a_n6_dual_a_n7_acc_match["acc_resn_x"] == a_n6_dual_a_n7_acc_match["don_resn_y"]) &
      (a_n6_dual_a_n7_acc_match["acc_resi_x"] == a_n6_dual_a_n7_acc_match["don_resi_y"]) &
      (a_n6_dual_a_n7_acc_match["acc_chain_x"] == a_n6_dual_a_n7_acc_match["don_chain_y"]) &
      ~(a_n6_dual_a_n7_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_dual_a_n7_acc_match["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_dual_a_n7_acc_match["acc_resn_x"] == a_n6_dual_a_n7_acc_match["don_resn_y"]) &
      (a_n6_dual_a_n7_acc_match["acc_resi_x"] == a_n6_dual_a_n7_acc_match["don_resi_y"]) &
      (a_n6_dual_a_n7_acc_match["acc_chain_x"] == a_n6_dual_a_n7_acc_match["don_chain_y"]) &
      ~(a_n6_dual_a_n7_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_dual_a_n7_acc_match["don_name_y"].isin(["O2'"])) &
      (a_n6_dual_a_n7_acc_match["acc_name_x"].isin(["O"])) &
      (a_n6_dual_a_n7_acc_match["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())

# Prepare dataframes with just residue information of single and dual H-bonding A(N6) that also accepts an H-bond at the
# N1, N3, or N7. Additionally, there is no overlap in the identity of the partner nucleobase, amino acid residue
# backbone, or amino acid residue side chain. In other words, the acceptor for the N6 H-bond donation and the donor for
# the N1, N3, or N7 H-bond acceptation belong to different partner entities.
a_n6_single_a_n1_acc_no_ol = pd.concat([a_n6_single_a_n1_acc_res, a_n6_single_a_n1_acc_ol]).drop_duplicates(keep=False)
a_n6_dual_a_n1_acc_no_ol = pd.concat([a_n6_dual_a_n1_acc_res, a_n6_dual_a_n1_acc_ol]).drop_duplicates(keep=False)
a_n6_single_a_n3_acc_no_ol = pd.concat([a_n6_single_a_n3_acc_res, a_n6_single_a_n3_acc_ol]).drop_duplicates(keep=False)
a_n6_dual_a_n3_acc_no_ol = pd.concat([a_n6_dual_a_n3_acc_res, a_n6_dual_a_n3_acc_ol]).drop_duplicates(keep=False)
a_n6_single_a_n7_acc_no_ol = pd.concat([a_n6_single_a_n7_acc_res, a_n6_single_a_n7_acc_ol]).drop_duplicates(keep=False)
a_n6_dual_a_n7_acc_no_ol = pd.concat([a_n6_dual_a_n7_acc_res, a_n6_dual_a_n7_acc_ol]).drop_duplicates(keep=False)

# Write data regarding A(N6) donors that also accept via the N1, N3, or N7 independent of same or different partner
# entity.
a_n6_donation = pd.DataFrame({
    "Atom": ["A(N1)", "A(N1)", "A(N1)", "A(N3)", "A(N3)", "A(N3)", "A(N7)", "A(N7)", "A(N7)"],
    "Type": ["No", "Single", "Dual", "No", "Single", "Dual", "No", "Single", "Dual"],
    "Occurrence": [a_n6_no_a_n1_acc_res.shape[0],
                   a_n6_single_a_n1_acc_res.shape[0],
                   a_n6_dual_a_n1_acc_res.shape[0],
                   a_n6_no_a_n3_acc_res.shape[0],
                   a_n6_single_a_n3_acc_res.shape[0],
                   a_n6_dual_a_n3_acc_res.shape[0],
                   a_n6_no_a_n7_acc_res.shape[0],
                   a_n6_single_a_n7_acc_res.shape[0],
                   a_n6_dual_a_n7_acc_res.shape[0]],
    "Total": [a_n6_no_res.shape[0],
              a_n6_single_res.shape[0],
              a_n6_dual_res.shape[0],
              a_n6_no_res.shape[0],
              a_n6_single_res.shape[0],
              a_n6_dual_res.shape[0],
              a_n6_no_res.shape[0],
              a_n6_single_res.shape[0],
              a_n6_dual_res.shape[0]]
})
a_n6_donation.to_csv("plots/frequency/a_n6_donation.csv", index=False)

# Write data regarding A(N6) donors that also accept via the N1, N3, or N7 with no overlap in partner entity.
a_n6_donation_no_ol = pd.DataFrame({
    "Atom": ["A(N1)", "A(N1)", "A(N1)", "A(N3)", "A(N3)", "A(N3)", "A(N7)", "A(N7)", "A(N7)"],
    "Type": ["No", "Single", "Dual", "No", "Single", "Dual", "No", "Single", "Dual"],
    "Occurrence": [a_n6_no_a_n1_acc_res.shape[0],
                   a_n6_single_a_n1_acc_no_ol.shape[0],
                   a_n6_dual_a_n1_acc_no_ol.shape[0],
                   a_n6_no_a_n3_acc_res.shape[0],
                   a_n6_single_a_n3_acc_no_ol.shape[0],
                   a_n6_dual_a_n3_acc_no_ol.shape[0],
                   a_n6_no_a_n7_acc_res.shape[0],
                   a_n6_single_a_n7_acc_no_ol.shape[0],
                   a_n6_dual_a_n7_acc_no_ol.shape[0]],
    "Total": [a_n6_no_res.shape[0],
              a_n6_single_res.shape[0],
              a_n6_dual_res.shape[0],
              a_n6_no_res.shape[0],
              a_n6_single_res.shape[0],
              a_n6_dual_res.shape[0],
              a_n6_no_res.shape[0],
              a_n6_single_res.shape[0],
              a_n6_dual_res.shape[0]]
})
a_n6_donation_no_ol.to_csv("plots/frequency/a_n6_donation_no_ol.csv", index=False)

# Prepare dataframes of A residues that do not donate via their N6 and that also accept via the N1, N3, or N7, including
# H-bond information related to the N1, N3, or N7 H-bond acceptation.
a_n6_no_a_n1_acc_hbond_to_n1 = (a_n6_no_a_n1_acc_res
                                .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]].eq(["N1", "A"])
                                       .all(axis='columns')], how='inner'))
a_n6_no_a_n3_acc_hbond_to_n3 = (a_n6_no_a_n3_acc_res
                                .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]].eq(["N3", "A"])
                                       .all(axis='columns')], how='inner'))
a_n6_no_a_n7_acc_hbond_to_n7 = (a_n6_no_a_n7_acc_res
                                .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]].eq(["N7", "A"])
                                       .all(axis='columns')], how='inner'))

# Prepare dataframes of A(N6) donors that also accept via the N1, N3, or N7 while only considering partner entities that
# do not overlap and which includes H-bond information related to the N1, N3, or N7 H-bond acceptation.
a_n6_single_a_n1_acc_no_ol_hbond_to_n1 = (a_n6_single_a_n1_acc_no_ol
                                          .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                           "chain": "acc_chain"})
                                          .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]]
                                                 .eq(["N1", "A"])
                                                 .all(axis='columns')], how='inner'))
a_n6_dual_a_n1_acc_no_ol_hbond_to_n1 = (a_n6_dual_a_n1_acc_no_ol
                                        .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                        .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]].eq(["N1", "A"])
                                               .all(axis='columns')], how='inner'))
a_n6_single_a_n3_acc_no_ol_hbond_to_n3 = (a_n6_single_a_n3_acc_no_ol
                                          .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                           "chain": "acc_chain"})
                                          .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]]
                                                 .eq(["N3", "A"])
                                                 .all(axis='columns')], how='inner'))
a_n6_dual_a_n3_acc_no_ol_hbond_to_n3 = (a_n6_dual_a_n3_acc_no_ol
                                        .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                        .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]].eq(["N3", "A"])
                                               .all(axis='columns')], how='inner'))
a_n6_single_a_n7_acc_no_ol_hbond_to_n7 = (a_n6_single_a_n7_acc_no_ol
                                          .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                           "chain": "acc_chain"})
                                          .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]]
                                                 .eq(["N7", "A"])
                                                 .all(axis='columns')], how='inner'))
a_n6_dual_a_n7_acc_no_ol_hbond_to_n7 = (a_n6_dual_a_n7_acc_no_ol
                                        .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                        .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]].eq(["N7", "A"])
                                               .all(axis='columns')], how='inner'))

# Write data on H-bond distance for donors to the N1, N3, or N7 of A residues categorized by the number of H-bonds they
# donate via their N6. Adenines that donate via their N6 have no overlap in partner entity when considering the donor
# to the N1, N3, or N7.
a_n6_donation_a_n1_acc_no_ol_hbond_to_n1 = pd.DataFrame({
    "Type": (["No"] * a_n6_no_a_n1_acc_hbond_to_n1.shape[0] +
             ["Single"] * a_n6_single_a_n1_acc_no_ol_hbond_to_n1.shape[0] +
             ["Dual"] * a_n6_dual_a_n1_acc_no_ol_hbond_to_n1.shape[0]),
    "Distance": (a_n6_no_a_n1_acc_hbond_to_n1["don_acc_distance"].to_list() +
                 a_n6_single_a_n1_acc_no_ol_hbond_to_n1["don_acc_distance"].to_list() +
                 a_n6_dual_a_n1_acc_no_ol_hbond_to_n1["don_acc_distance"].to_list())
})
a_n6_donation_a_n1_acc_no_ol_hbond_to_n1.to_csv("plots/dist/a_n1_acc/a_n6_donation_a_n1_acc_no_ol_hbond_to_n1.csv",
                                                index=False)
a_n6_donation_a_n3_acc_no_ol_hbond_to_n3 = pd.DataFrame({
    "Type": (["No"] * a_n6_no_a_n3_acc_hbond_to_n3.shape[0] +
             ["Single"] * a_n6_single_a_n3_acc_no_ol_hbond_to_n3.shape[0] +
             ["Dual"] * a_n6_dual_a_n3_acc_no_ol_hbond_to_n3.shape[0]),
    "Distance": (a_n6_no_a_n3_acc_hbond_to_n3["don_acc_distance"].to_list() +
                 a_n6_single_a_n3_acc_no_ol_hbond_to_n3["don_acc_distance"].to_list() +
                 a_n6_dual_a_n3_acc_no_ol_hbond_to_n3["don_acc_distance"].to_list())
})
a_n6_donation_a_n3_acc_no_ol_hbond_to_n3.to_csv("plots/dist/a_n3_acc/a_n6_donation_a_n3_acc_no_ol_hbond_to_n3.csv",
                                                index=False)
a_n6_donation_a_n7_acc_no_ol_hbond_to_n7 = pd.DataFrame({
    "Type": (["No"] * a_n6_no_a_n7_acc_hbond_to_n7.shape[0] +
             ["Single"] * a_n6_single_a_n7_acc_no_ol_hbond_to_n7.shape[0] +
             ["Dual"] * a_n6_dual_a_n7_acc_no_ol_hbond_to_n7.shape[0]),
    "Distance": (a_n6_no_a_n7_acc_hbond_to_n7["don_acc_distance"].to_list() +
                 a_n6_single_a_n7_acc_no_ol_hbond_to_n7["don_acc_distance"].to_list() +
                 a_n6_dual_a_n7_acc_no_ol_hbond_to_n7["don_acc_distance"].to_list())
})
a_n6_donation_a_n7_acc_no_ol_hbond_to_n7.to_csv("plots/dist/a_n7_acc/a_n6_donation_a_n7_acc_no_ol_hbond_to_n7.csv",
                                                index=False)

# Prepare dataframes of single and dual H-bonding A(N6) residues that only donate to OD1 or OD2 of Asp, OE1 or OE2 of
# Glu, or OP1 or OP2 of nucleic acids and that accept an H-bond at the N1. The dataframes include H-bond information
# related to the N6 H-bond donation.
find_neg_single_n1 = a_n6_single_a_n1_acc_hbond_from_n6
find_neg_single_n1["neg_acc"] = pd.Series((a_n6_single_a_n1_acc_hbond_from_n6["acc_resn"].isin(neg_acc_resn)) &
                                          (a_n6_single_a_n1_acc_hbond_from_n6["acc_name"].isin(neg_acc_name)))
a_n6_single_a_n1_acc_hbond_from_n6_neg = (find_neg_single_n1[find_neg_single_n1
                                          .groupby(["don_index", "eq_class_members"])["neg_acc"]
                                          .transform("all")].drop(columns=["neg_acc"]))
find_neg_dual = a_n6_dual_a_n1_acc_hbond_from_n6
find_neg_dual["neg_acc"] = pd.Series((a_n6_dual_a_n1_acc_hbond_from_n6["acc_resn"].isin(neg_acc_resn)) &
                                     (a_n6_dual_a_n1_acc_hbond_from_n6["acc_name"].isin(neg_acc_name)))
a_n6_dual_a_n1_acc_hbond_from_n6_neg = (find_neg_dual[find_neg_dual.groupby(["don_index", "eq_class_members"])["neg_acc"]
                                        .transform("all")].drop(columns=["neg_acc"]))

# Prepare dataframes of single and dual H-bonding A(N6) residues that only donate to OD1 or OD2 of Asp, OE1 or OE2 of
# Glu, or OP1 or OP2 of nucleic acids and that accept an H-bond at the N3. The dataframes include H-bond information
# related to the N6 H-bond donation.
find_neg_single_n3 = a_n6_single_a_n3_acc_hbond_from_n6
find_neg_single_n3["neg_acc"] = pd.Series((a_n6_single_a_n3_acc_hbond_from_n6["acc_resn"].isin(neg_acc_resn)) &
                                          (a_n6_single_a_n3_acc_hbond_from_n6["acc_name"].isin(neg_acc_name)))
a_n6_single_a_n3_acc_hbond_from_n6_neg = (find_neg_single_n3[find_neg_single_n3
                                          .groupby(["don_index", "eq_class_members"])["neg_acc"]
                                          .transform("all")].drop(columns=["neg_acc"]))
find_neg_dual = a_n6_dual_a_n3_acc_hbond_from_n6
find_neg_dual["neg_acc"] = pd.Series((a_n6_dual_a_n3_acc_hbond_from_n6["acc_resn"].isin(neg_acc_resn)) &
                                     (a_n6_dual_a_n3_acc_hbond_from_n6["acc_name"].isin(neg_acc_name)))
a_n6_dual_a_n3_acc_hbond_from_n6_neg = (find_neg_dual[find_neg_dual.groupby(["don_index", "eq_class_members"])["neg_acc"]
                                        .transform("all")].drop(columns=["neg_acc"]))

# Prepare dataframes of single and dual H-bonding A(N6) residues that only donate to OD1 or OD2 of Asp, OE1 or OE2 of
# Glu, or OP1 or OP2 of nucleic acids and that accept an H-bond at the N7. The dataframes include H-bond information
# related to the N6 H-bond donation.
find_neg_single_n7 = a_n6_single_a_n7_acc_hbond_from_n6
find_neg_single_n7["neg_acc"] = pd.Series((a_n6_single_a_n7_acc_hbond_from_n6["acc_resn"].isin(neg_acc_resn)) &
                                          (a_n6_single_a_n7_acc_hbond_from_n6["acc_name"].isin(neg_acc_name)))
a_n6_single_a_n7_acc_hbond_from_n6_neg = (find_neg_single_n7[find_neg_single_n7
                                          .groupby(["don_index", "eq_class_members"])["neg_acc"]
                                          .transform("all")].drop(columns=["neg_acc"]))
find_neg_dual = a_n6_dual_a_n7_acc_hbond_from_n6
find_neg_dual["neg_acc"] = pd.Series((a_n6_dual_a_n7_acc_hbond_from_n6["acc_resn"].isin(neg_acc_resn)) &
                                     (a_n6_dual_a_n7_acc_hbond_from_n6["acc_name"].isin(neg_acc_name)))
a_n6_dual_a_n7_acc_hbond_from_n6_neg = (find_neg_dual[find_neg_dual.groupby(["don_index", "eq_class_members"])["neg_acc"]
                                        .transform("all")].drop(columns=["neg_acc"]))

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at the N1 with just residue information. The N6
# only donates to OD1 or OD2 of Asp, OE1 or OE2 of Glu, or OP1 or OP2 of nucleic acids.
a_n6_single_a_n1_acc_res_neg = ((a_n6_single_a_n1_acc_hbond_from_n6_neg
                                .loc[:, ["eq_class_members", "don_resn", "don_resi", "don_chain"]]
                                .drop_duplicates())
                                .rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"}))
a_n6_dual_a_n1_acc_res_neg = ((a_n6_dual_a_n1_acc_hbond_from_n6_neg
                              .loc[:, ["eq_class_members", "don_resn", "don_resi", "don_chain"]]
                              .drop_duplicates())
                              .rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"}))

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at the N3 with just residue information. The N6
# only donates to OD1 or OD2 of Asp, OE1 or OE2 of Glu, or OP1 or OP2 of nucleic acids.
a_n6_single_a_n3_acc_res_neg = ((a_n6_single_a_n3_acc_hbond_from_n6_neg
                                .loc[:, ["eq_class_members", "don_resn", "don_resi", "don_chain"]]
                                .drop_duplicates())
                                .rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"}))
a_n6_dual_a_n3_acc_res_neg = ((a_n6_dual_a_n3_acc_hbond_from_n6_neg
                              .loc[:, ["eq_class_members", "don_resn", "don_resi", "don_chain"]]
                              .drop_duplicates())
                              .rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"}))

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at the N7 with just residue information. The N6
# only donates to OD1 or OD2 of Asp, OE1 or OE2 of Glu, or OP1 or OP2 of nucleic acids.
a_n6_single_a_n7_acc_res_neg = ((a_n6_single_a_n7_acc_hbond_from_n6_neg
                                .loc[:, ["eq_class_members", "don_resn", "don_resi", "don_chain"]]
                                .drop_duplicates())
                                .rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"}))
a_n6_dual_a_n7_acc_res_neg = ((a_n6_dual_a_n7_acc_hbond_from_n6_neg
                              .loc[:, ["eq_class_members", "don_resn", "don_resi", "don_chain"]]
                              .drop_duplicates())
                              .rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"}))

# Prepare dataframes of single and dual H-bonding A(N6) residues that only donate to OD1 or OD2 of Asp, OE1 or OE2 of
# Glu, or OP1 or OP2 of nucleic acids and that accept an H-bond at the N1. The dataframes include H-bond information
# related to the N6 H-bond donation and the N1 H-bond acceptation.
a_n6_single_a_n1_acc_match_neg = (a_n6_single_a_n1_acc_hbond_from_n6_neg
                                  .merge(a_n6_single_a_n1_acc_hbond_to_n1,
                                         left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                         right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"], how='inner'))
a_n6_dual_a_n1_acc_match_neg = (a_n6_dual_a_n1_acc_hbond_from_n6_neg
                                .merge(a_n6_dual_a_n1_acc_hbond_to_n1,
                                       left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                       right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"], how='inner'))

# Prepare dataframes of single and dual H-bonding A(N6) residues that only donate to OD1 or OD2 of Asp, OE1 or OE2 of
# Glu, or OP1 or OP2 of nucleic acids and that accept an H-bond at the N3. The dataframes include H-bond information
# related to the N6 H-bond donation and the N3 H-bond acceptation.
a_n6_single_a_n3_acc_match_neg = (a_n6_single_a_n3_acc_hbond_from_n6_neg
                                  .merge(a_n6_single_a_n3_acc_hbond_to_n3,
                                         left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                         right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"], how='inner'))
a_n6_dual_a_n3_acc_match_neg = (a_n6_dual_a_n3_acc_hbond_from_n6_neg
                                .merge(a_n6_dual_a_n3_acc_hbond_to_n3,
                                       left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                       right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"], how='inner'))

# Prepare dataframes of single and dual H-bonding A(N6) residues that only donate to OD1 or OD2 of Asp, OE1 or OE2 of
# Glu, or OP1 or OP2 of nucleic acids and that accept an H-bond at the N7. The dataframes include H-bond information
# related to the N6 H-bond donation and the N7 H-bond acceptation.
a_n6_single_a_n7_acc_match_neg = (a_n6_single_a_n7_acc_hbond_from_n6_neg
                                  .merge(a_n6_single_a_n7_acc_hbond_to_n7,
                                         left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                         right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"], how='inner'))
a_n6_dual_a_n7_acc_match_neg = (a_n6_dual_a_n7_acc_hbond_from_n6_neg
                                .merge(a_n6_dual_a_n7_acc_hbond_to_n7,
                                       left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                       right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"], how='inner'))

# Prepare dataframes with just residue information of single and dual H-bonding A(N6) that also accepts an H-bond at the
# N1, N3, or N7. The N6 only donates to OD1 or OD2 of Asp, OE1 or OE2 of Glu, or OP1 or OP2 of nucleic acids.
# Additionally, there is overlap in the identity of the partner nucleobase, amino acid residue backbone, or amino acid
# residue side chain, where the N6 donates an H-bond and the N1, N3, or N7 accepts an H-bond from the same partner
# entity.
a_n6_single_a_n1_acc_ol_neg = (a_n6_single_a_n1_acc_match_neg[
    (((a_n6_single_a_n1_acc_match_neg["acc_resn_x"] == a_n6_single_a_n1_acc_match_neg["don_resn_y"]) &
      (a_n6_single_a_n1_acc_match_neg["acc_resi_x"] == a_n6_single_a_n1_acc_match_neg["don_resi_y"]) &
      (a_n6_single_a_n1_acc_match_neg["acc_chain_x"] == a_n6_single_a_n1_acc_match_neg["don_chain_y"]) &
      ~(a_n6_single_a_n1_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_single_a_n1_acc_match_neg["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_single_a_n1_acc_match_neg["acc_resn_x"] == a_n6_single_a_n1_acc_match_neg["don_resn_y"]) &
      (a_n6_single_a_n1_acc_match_neg["acc_resi_x"] == a_n6_single_a_n1_acc_match_neg["don_resi_y"]) &
      (a_n6_single_a_n1_acc_match_neg["acc_chain_x"] == a_n6_single_a_n1_acc_match_neg["don_chain_y"]) &
      ~(a_n6_single_a_n1_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_single_a_n1_acc_match_neg["don_name_y"].isin(["O2'"])) &
      (a_n6_single_a_n1_acc_match_neg["acc_name_x"].isin(["O"])) &
      (a_n6_single_a_n1_acc_match_neg["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())
a_n6_dual_a_n1_acc_ol_neg = (a_n6_dual_a_n1_acc_match_neg[
    (((a_n6_dual_a_n1_acc_match_neg["acc_resn_x"] == a_n6_dual_a_n1_acc_match_neg["don_resn_y"]) &
      (a_n6_dual_a_n1_acc_match_neg["acc_resi_x"] == a_n6_dual_a_n1_acc_match_neg["don_resi_y"]) &
      (a_n6_dual_a_n1_acc_match_neg["acc_chain_x"] == a_n6_dual_a_n1_acc_match_neg["don_chain_y"]) &
      ~(a_n6_dual_a_n1_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_dual_a_n1_acc_match_neg["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_dual_a_n1_acc_match_neg["acc_resn_x"] == a_n6_dual_a_n1_acc_match_neg["don_resn_y"]) &
      (a_n6_dual_a_n1_acc_match_neg["acc_resi_x"] == a_n6_dual_a_n1_acc_match_neg["don_resi_y"]) &
      (a_n6_dual_a_n1_acc_match_neg["acc_chain_x"] == a_n6_dual_a_n1_acc_match_neg["don_chain_y"]) &
      ~(a_n6_dual_a_n1_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_dual_a_n1_acc_match_neg["don_name_y"].isin(["O2'"])) &
      (a_n6_dual_a_n1_acc_match_neg["acc_name_x"].isin(["O"])) &
      (a_n6_dual_a_n1_acc_match_neg["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())
a_n6_single_a_n3_acc_ol_neg = (a_n6_single_a_n3_acc_match_neg[
    (((a_n6_single_a_n3_acc_match_neg["acc_resn_x"] == a_n6_single_a_n3_acc_match_neg["don_resn_y"]) &
      (a_n6_single_a_n3_acc_match_neg["acc_resi_x"] == a_n6_single_a_n3_acc_match_neg["don_resi_y"]) &
      (a_n6_single_a_n3_acc_match_neg["acc_chain_x"] == a_n6_single_a_n3_acc_match_neg["don_chain_y"]) &
      ~(a_n6_single_a_n3_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_single_a_n3_acc_match_neg["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_single_a_n3_acc_match_neg["acc_resn_x"] == a_n6_single_a_n3_acc_match_neg["don_resn_y"]) &
      (a_n6_single_a_n3_acc_match_neg["acc_resi_x"] == a_n6_single_a_n3_acc_match_neg["don_resi_y"]) &
      (a_n6_single_a_n3_acc_match_neg["acc_chain_x"] == a_n6_single_a_n3_acc_match_neg["don_chain_y"]) &
      ~(a_n6_single_a_n3_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_single_a_n3_acc_match_neg["don_name_y"].isin(["O2'"])) &
      (a_n6_single_a_n3_acc_match_neg["acc_name_x"].isin(["O"])) &
      (a_n6_single_a_n3_acc_match_neg["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())
a_n6_dual_a_n3_acc_ol_neg = (a_n6_dual_a_n3_acc_match_neg[
    (((a_n6_dual_a_n3_acc_match_neg["acc_resn_x"] == a_n6_dual_a_n3_acc_match_neg["don_resn_y"]) &
      (a_n6_dual_a_n3_acc_match_neg["acc_resi_x"] == a_n6_dual_a_n3_acc_match_neg["don_resi_y"]) &
      (a_n6_dual_a_n3_acc_match_neg["acc_chain_x"] == a_n6_dual_a_n3_acc_match_neg["don_chain_y"]) &
      ~(a_n6_dual_a_n3_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_dual_a_n3_acc_match_neg["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_dual_a_n3_acc_match_neg["acc_resn_x"] == a_n6_dual_a_n3_acc_match_neg["don_resn_y"]) &
      (a_n6_dual_a_n3_acc_match_neg["acc_resi_x"] == a_n6_dual_a_n3_acc_match_neg["don_resi_y"]) &
      (a_n6_dual_a_n3_acc_match_neg["acc_chain_x"] == a_n6_dual_a_n3_acc_match_neg["don_chain_y"]) &
      ~(a_n6_dual_a_n3_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_dual_a_n3_acc_match_neg["don_name_y"].isin(["O2'"])) &
      (a_n6_dual_a_n3_acc_match_neg["acc_name_x"].isin(["O"])) &
      (a_n6_dual_a_n3_acc_match_neg["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())
a_n6_single_a_n7_acc_ol_neg = (a_n6_single_a_n7_acc_match_neg[
    (((a_n6_single_a_n7_acc_match_neg["acc_resn_x"] == a_n6_single_a_n7_acc_match_neg["don_resn_y"]) &
      (a_n6_single_a_n7_acc_match_neg["acc_resi_x"] == a_n6_single_a_n7_acc_match_neg["don_resi_y"]) &
      (a_n6_single_a_n7_acc_match_neg["acc_chain_x"] == a_n6_single_a_n7_acc_match_neg["don_chain_y"]) &
      ~(a_n6_single_a_n7_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_single_a_n7_acc_match_neg["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_single_a_n7_acc_match_neg["acc_resn_x"] == a_n6_single_a_n7_acc_match_neg["don_resn_y"]) &
      (a_n6_single_a_n7_acc_match_neg["acc_resi_x"] == a_n6_single_a_n7_acc_match_neg["don_resi_y"]) &
      (a_n6_single_a_n7_acc_match_neg["acc_chain_x"] == a_n6_single_a_n7_acc_match_neg["don_chain_y"]) &
      ~(a_n6_single_a_n7_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_single_a_n7_acc_match_neg["don_name_y"].isin(["O2'"])) &
      (a_n6_single_a_n7_acc_match_neg["acc_name_x"].isin(["O"])) &
      (a_n6_single_a_n7_acc_match_neg["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())
a_n6_dual_a_n7_acc_ol_neg = (a_n6_dual_a_n7_acc_match_neg[
    (((a_n6_dual_a_n7_acc_match_neg["acc_resn_x"] == a_n6_dual_a_n7_acc_match_neg["don_resn_y"]) &
      (a_n6_dual_a_n7_acc_match_neg["acc_resi_x"] == a_n6_dual_a_n7_acc_match_neg["don_resi_y"]) &
      (a_n6_dual_a_n7_acc_match_neg["acc_chain_x"] == a_n6_dual_a_n7_acc_match_neg["don_chain_y"]) &
      ~(a_n6_dual_a_n7_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_dual_a_n7_acc_match_neg["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_dual_a_n7_acc_match_neg["acc_resn_x"] == a_n6_dual_a_n7_acc_match_neg["don_resn_y"]) &
      (a_n6_dual_a_n7_acc_match_neg["acc_resi_x"] == a_n6_dual_a_n7_acc_match_neg["don_resi_y"]) &
      (a_n6_dual_a_n7_acc_match_neg["acc_chain_x"] == a_n6_dual_a_n7_acc_match_neg["don_chain_y"]) &
      ~(a_n6_dual_a_n7_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_dual_a_n7_acc_match_neg["don_name_y"].isin(["O2'"])) &
      (a_n6_dual_a_n7_acc_match_neg["acc_name_x"].isin(["O"])) &
      (a_n6_dual_a_n7_acc_match_neg["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())

# Prepare dataframes with just residue information of single and dual H-bonding A(N6) that also accepts an H-bond at the
# N1, N3, or N7. The N6 only donates to OD1 or OD2 of Asp, OE1 or OE2 of Glu, or OP1 or OP2 of nucleic acids.
# Additionally, there is no overlap in the identity of the partner nucleobase, amino acid residue backbone, or amino
# acid residue side chain. In other words, the acceptor for the N6 H-bond donation and the donor for the N1, N3, or N7
# H-bond acceptation belong to different partner entities.
a_n6_single_a_n1_acc_no_ol_neg = (pd.concat([a_n6_single_a_n1_acc_res_neg, a_n6_single_a_n1_acc_ol_neg])
                                  .drop_duplicates(keep=False))
a_n6_dual_a_n1_acc_no_ol_neg = (pd.concat([a_n6_dual_a_n1_acc_res_neg, a_n6_dual_a_n1_acc_ol_neg])
                                .drop_duplicates(keep=False))
a_n6_single_a_n3_acc_no_ol_neg = (pd.concat([a_n6_single_a_n3_acc_res_neg, a_n6_single_a_n3_acc_ol_neg])
                                  .drop_duplicates(keep=False))
a_n6_dual_a_n3_acc_no_ol_neg = (pd.concat([a_n6_dual_a_n3_acc_res_neg, a_n6_dual_a_n3_acc_ol_neg])
                                .drop_duplicates(keep=False))
a_n6_single_a_n7_acc_no_ol_neg = (pd.concat([a_n6_single_a_n7_acc_res_neg, a_n6_single_a_n7_acc_ol_neg])
                                  .drop_duplicates(keep=False))
a_n6_dual_a_n7_acc_no_ol_neg = (pd.concat([a_n6_dual_a_n7_acc_res_neg, a_n6_dual_a_n7_acc_ol_neg])
                                .drop_duplicates(keep=False))

# Prepare dataframes of A(N6) donors that also accept via the N1, N3, or N7 while only considering partner entities that
# do not overlap and which includes H-bond information related to the N1, N3, or N7 H-bond acceptation. The N6 only
# donates to OD1 or OD2 of Asp, OE1 or OE2 of Glu, or OP1 or OP2 of nucleic acids.
a_n6_single_a_n1_acc_no_ol_hbond_to_n1_neg = (a_n6_single_a_n1_acc_no_ol_neg
                                              .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                               "chain": "acc_chain"})
                                              .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]]
                                                     .eq(["N1", "A"])
                                                     .all(axis='columns')], how='inner'))
a_n6_dual_a_n1_acc_no_ol_hbond_to_n1_neg = (a_n6_dual_a_n1_acc_no_ol_neg
                                            .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                             "chain": "acc_chain"})
                                            .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]]
                                                   .eq(["N1", "A"])
                                                   .all(axis='columns')], how='inner'))
a_n6_single_a_n3_acc_no_ol_hbond_to_n3_neg = (a_n6_single_a_n3_acc_no_ol_neg
                                              .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                               "chain": "acc_chain"})
                                              .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]]
                                                     .eq(["N3", "A"])
                                                     .all(axis='columns')], how='inner'))
a_n6_dual_a_n3_acc_no_ol_hbond_to_n3_neg = (a_n6_dual_a_n3_acc_no_ol_neg
                                            .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                             "chain": "acc_chain"})
                                            .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]]
                                                   .eq(["N3", "A"])
                                                   .all(axis='columns')], how='inner'))
a_n6_single_a_n7_acc_no_ol_hbond_to_n7_neg = (a_n6_single_a_n7_acc_no_ol_neg
                                              .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                               "chain": "acc_chain"})
                                              .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]]
                                                     .eq(["N7", "A"])
                                                     .all(axis='columns')], how='inner'))
a_n6_dual_a_n7_acc_no_ol_hbond_to_n7_neg = (a_n6_dual_a_n7_acc_no_ol_neg
                                            .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                             "chain": "acc_chain"})
                                            .merge(acc_hbonds_nr[acc_hbonds_nr[["acc_name", "acc_resn"]]
                                                   .eq(["N7", "A"])
                                                   .all(axis='columns')], how='inner'))

# Write data on H-bond distance for donors to the N1, N3, or N7 of A residues that also donate via their
# N6 with no overlap in partner entity. Include data where the N6 donates to all acceptor types and where it only
# donates to OD1 or OD2 of Asp, OE1 or OE2 of Glu, or OP1 or OP2 of nucleic acids.
a_n6_donation_a_n1_acc_no_ol_hbond_to_n1_all_neg = pd.DataFrame({
    "Acceptor": (["All"] * a_n6_single_a_n1_acc_no_ol_hbond_to_n1.shape[0] +
                 ["All"] * a_n6_dual_a_n1_acc_no_ol_hbond_to_n1.shape[0] +
                 ["Negative"] * a_n6_single_a_n1_acc_no_ol_hbond_to_n1_neg.shape[0] +
                 ["Negative"] * a_n6_dual_a_n1_acc_no_ol_hbond_to_n1_neg.shape[0]),
    "Type": (["Single"] * a_n6_single_a_n1_acc_no_ol_hbond_to_n1.shape[0] +
             ["Dual"] * a_n6_dual_a_n1_acc_no_ol_hbond_to_n1.shape[0] +
             ["Single"] * a_n6_single_a_n1_acc_no_ol_hbond_to_n1_neg.shape[0] +
             ["Dual"] * a_n6_dual_a_n1_acc_no_ol_hbond_to_n1_neg.shape[0]),
    "Distance": (a_n6_single_a_n1_acc_no_ol_hbond_to_n1["don_acc_distance"].to_list() +
                 a_n6_dual_a_n1_acc_no_ol_hbond_to_n1["don_acc_distance"].to_list() +
                 a_n6_single_a_n1_acc_no_ol_hbond_to_n1_neg["don_acc_distance"].to_list() +
                 a_n6_dual_a_n1_acc_no_ol_hbond_to_n1_neg["don_acc_distance"].to_list())
})
a_n6_donation_a_n1_acc_no_ol_hbond_to_n1_all_neg.to_csv(
    "plots/dist/a_n1_acc/a_n6_donation_a_n1_acc_no_ol_hbond_to_n1_all_neg.csv", index=False)
a_n6_donation_a_n3_acc_no_ol_hbond_to_n3_all_neg = pd.DataFrame({
    "Acceptor": (["All"] * a_n6_single_a_n3_acc_no_ol_hbond_to_n3.shape[0] +
                 ["All"] * a_n6_dual_a_n3_acc_no_ol_hbond_to_n3.shape[0] +
                 ["Negative"] * a_n6_single_a_n3_acc_no_ol_hbond_to_n3_neg.shape[0] +
                 ["Negative"] * a_n6_dual_a_n3_acc_no_ol_hbond_to_n3_neg.shape[0]),
    "Type": (["Single"] * a_n6_single_a_n3_acc_no_ol_hbond_to_n3.shape[0] +
             ["Dual"] * a_n6_dual_a_n3_acc_no_ol_hbond_to_n3.shape[0] +
             ["Single"] * a_n6_single_a_n3_acc_no_ol_hbond_to_n3_neg.shape[0] +
             ["Dual"] * a_n6_dual_a_n3_acc_no_ol_hbond_to_n3_neg.shape[0]),
    "Distance": (a_n6_single_a_n3_acc_no_ol_hbond_to_n3["don_acc_distance"].to_list() +
                 a_n6_dual_a_n3_acc_no_ol_hbond_to_n3["don_acc_distance"].to_list() +
                 a_n6_single_a_n3_acc_no_ol_hbond_to_n3_neg["don_acc_distance"].to_list() +
                 a_n6_dual_a_n3_acc_no_ol_hbond_to_n3_neg["don_acc_distance"].to_list())
})
a_n6_donation_a_n3_acc_no_ol_hbond_to_n3_all_neg.to_csv(
    "plots/dist/a_n3_acc/a_n6_donation_a_n3_acc_no_ol_hbond_to_n3_all_neg.csv", index=False)
a_n6_donation_a_n7_acc_no_ol_hbond_to_n7_all_neg = pd.DataFrame({
    "Acceptor": (["All"] * a_n6_single_a_n7_acc_no_ol_hbond_to_n7.shape[0] +
                 ["All"] * a_n6_dual_a_n7_acc_no_ol_hbond_to_n7.shape[0] +
                 ["Negative"] * a_n6_single_a_n7_acc_no_ol_hbond_to_n7_neg.shape[0] +
                 ["Negative"] * a_n6_dual_a_n7_acc_no_ol_hbond_to_n7_neg.shape[0]),
    "Type": (["Single"] * a_n6_single_a_n7_acc_no_ol_hbond_to_n7.shape[0] +
             ["Dual"] * a_n6_dual_a_n7_acc_no_ol_hbond_to_n7.shape[0] +
             ["Single"] * a_n6_single_a_n7_acc_no_ol_hbond_to_n7_neg.shape[0] +
             ["Dual"] * a_n6_dual_a_n7_acc_no_ol_hbond_to_n7_neg.shape[0]),
    "Distance": (a_n6_single_a_n7_acc_no_ol_hbond_to_n7["don_acc_distance"].to_list() +
                 a_n6_dual_a_n7_acc_no_ol_hbond_to_n7["don_acc_distance"].to_list() +
                 a_n6_single_a_n7_acc_no_ol_hbond_to_n7_neg["don_acc_distance"].to_list() +
                 a_n6_dual_a_n7_acc_no_ol_hbond_to_n7_neg["don_acc_distance"].to_list())
})
a_n6_donation_a_n7_acc_no_ol_hbond_to_n7_all_neg.to_csv(
    "plots/dist/a_n7_acc/a_n6_donation_a_n7_acc_no_ol_hbond_to_n7_all_neg.csv", index=False)

# TODO consider revising collect_data.py to only consider nucleobases with all heavy atoms, see resn A and resi 2150 and chain A eq_class_members NR_3.0_90633.21 for an example of an incomplete nucleobase, leading to the following evaluation to be False
# print((a_n6_no_a_n7_acc_res["eq_class_members"].size + a_n6_single_a_n7_acc_res["eq_class_members"].size + a_n6_dual_a_n7_acc_res["eq_class_members"].size) == a_n7_acc_res["eq_class_members"].size)
