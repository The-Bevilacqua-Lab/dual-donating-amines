"""
This script will eventually do something.
"""

import os
import pandas as pd
import numpy as np
import scipy

# the working directory should house the combined/ folder
working_dir = os.getcwd()

# create the combined dataframes
don_hbonds_nr_c = pd.read_csv(working_dir + "/combined/don_hbonds_nr_c.csv",
                              dtype={"don_resi": "object", "acc_resi": "object"})
prot_don_hbonds_nr_c = pd.read_csv(working_dir + "/combined/prot_don_hbonds_nr_c.csv",
                                   dtype={"don_resi": "object", "acc_resi": "object"})
acc_hbonds_nr_c = pd.read_csv(working_dir + "/combined/acc_hbonds_nr_c.csv",
                              dtype={"don_resi": "object", "acc_resi": "object"})
deprot_acc_hbonds_nr_c = pd.read_csv(working_dir + "/combined/deprot_acc_hbonds_nr_c.csv",
                                     dtype={"don_resi": "object", "acc_resi": "object"})
nuc_data_c = pd.read_csv(working_dir + "/combined/nuc_data_c.csv")
single_don_exo_amines_c = pd.read_csv(working_dir + "/combined/single_don_exo_amines_c.csv")
dual_don_exo_amines_c = pd.read_csv(working_dir + "/combined/dual_don_exo_amines_c.csv")

# prepare dataframes of exocyclic amines for each A, C, and G residue
a_exo_amines = nuc_data_c[nuc_data_c["atom_name"] == "N6"]

# prepare dataframes of single and dual H-bonding A(N6) with H-bond information
a_n6_single_hbond = single_don_exo_amines_c.merge(don_hbonds_nr_c[don_hbonds_nr_c["don_resn"] == "A"],
                                                on=["don_index", "eq_class"], how='inner')
a_n6_dual_hbond = dual_don_exo_amines_c.merge(don_hbonds_nr_c[don_hbonds_nr_c["don_resn"] == "A"],
                                            on=["don_index", "eq_class"], how='inner')

# prepare dataframes of no, single, and dual H-bonding A(N6) with just residue information
a_n6_single_res = (a_n6_single_hbond.loc[:, ["eq_class", "don_resn", "don_resi", "don_chain"]]
                   .drop_duplicates()).rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"})
a_n6_dual_res = (a_n6_dual_hbond.loc[:, ["eq_class", "don_resn", "don_resi", "don_chain"]]
                 .drop_duplicates()).rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"})
a_n6_no_res = (pd.concat([a_n6_single_res,
                         a_n6_dual_res,
                         a_exo_amines.loc[:, ["eq_class", "resn", "resi", "chain"]]])
               .drop_duplicates(keep=False))

# TODO re-evaluate the following two variables at some point if I am going to keep them
# prepare a list of single and dual H-bonding A residues that are involved in a canonical AU base pair
au_bp_single = (a_n6_single_hbond[a_n6_single_hbond[["acc_name", "acc_resn"]].eq(["O4", "U"]).all(axis='columns')]
                .merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["don_name", "don_resn", "acc_name"]].eq(["N3", "U", "N1"])
                       .all(axis='columns')],
                       left_on=["don_resn", "don_resi", "don_chain", "acc_resn", "acc_resi", "acc_chain", "eq_class"],
                       right_on=["acc_resn", "acc_resi", "acc_chain", "don_resn", "don_resi", "don_chain", "eq_class"],
                       how='inner'))
au_bp_dual = (a_n6_dual_hbond[a_n6_dual_hbond[["acc_name", "acc_resn"]].eq(["O4", "U"]).all(axis='columns')]
              .merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["don_name", "don_resn", "acc_name"]].eq(["N3", "U", "N1"])
                     .all(axis='columns')],
                     left_on=["don_resn", "don_resi", "don_chain", "acc_resn", "acc_resi", "acc_chain", "eq_class"],
                     right_on=["acc_resn", "acc_resi", "acc_chain", "don_resn", "don_resi", "don_chain", "eq_class"],
                     how='inner'))

# prepare dataframes with just residue information of A residues that accept H-bonds at their N1, N3, or N7
a_n1_acc_res = (acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]].eq(["N1", "A"]).all(axis='columns')]
                .loc[:, ["eq_class", "acc_resn", "acc_resi", "acc_chain"]].drop_duplicates()
                .rename(columns={"acc_resn": "resn", "acc_resi": "resi", "acc_chain": "chain"}))
a_n3_acc_res = (acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]].eq(["N3", "A"]).all(axis='columns')]
                .loc[:, ["eq_class", "acc_resn", "acc_resi", "acc_chain"]].drop_duplicates()
                .rename(columns={"acc_resn": "resn", "acc_resi": "resi", "acc_chain": "chain"}))
a_n7_acc_res = (acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]].eq(["N7", "A"]).all(axis='columns')]
                .loc[:, ["eq_class", "acc_resn", "acc_resi", "acc_chain"]].drop_duplicates()
                .rename(columns={"acc_resn": "resn", "acc_resi": "resi", "acc_chain": "chain"}))

# prepare dataframes of no, single, and dual H-bonding A(N6) that also accept at the N1 with just residue information
a_n6_no_a_n1_acc_res = a_n6_no_res.merge(a_n1_acc_res, how='inner')
a_n6_single_a_n1_acc_res = a_n6_single_res.merge(a_n1_acc_res, how='inner')
a_n6_dual_a_n1_acc_res = a_n6_dual_res.merge(a_n1_acc_res, how='inner')

# prepare dataframes of no, single, and dual H-bonding A(N6) that also accept at the N3 with just residue information
a_n6_no_a_n3_acc_res = a_n6_no_res.merge(a_n3_acc_res, how='inner')
a_n6_single_a_n3_acc_res = a_n6_single_res.merge(a_n3_acc_res, how='inner')
a_n6_dual_a_n3_acc_res = a_n6_dual_res.merge(a_n3_acc_res, how='inner')

# prepare dataframes of no, single, and dual H-bonding A(N6) that also accept at the N7 with just residue information
a_n6_no_a_n7_acc_res = a_n6_no_res.merge(a_n7_acc_res, how='inner')
a_n6_single_a_n7_acc_res = a_n6_single_res.merge(a_n7_acc_res, how='inner')
a_n6_dual_a_n7_acc_res = a_n6_dual_res.merge(a_n7_acc_res, how='inner')

# prepare dataframes of single and dual H-bonding A(N6) that also accept at N1 and which include H-bond information
# related to the N6 H-bond donation
a_n6_single_a_n1_acc_hbond_from_n6 = (a_n6_single_a_n1_acc_res
                                      .rename(columns={"resn": "don_resn", "resi": "don_resi", "chain": "don_chain"})
                                      .merge(a_n6_single_hbond, how='inner'))
a_n6_dual_a_n1_acc_hbond_from_n6 = (a_n6_dual_a_n1_acc_res
                                    .rename(columns={"resn": "don_resn", "resi": "don_resi", "chain": "don_chain"})
                                    .merge(a_n6_dual_hbond, how='inner'))

# prepare dataframes of single and dual H-bonding A(N6) that also accept at N3 and which include H-bond information
# related to the N6 H-bond donation
a_n6_single_a_n3_acc_hbond_from_n6 = (a_n6_single_a_n3_acc_res
                                      .rename(columns={"resn": "don_resn", "resi": "don_resi", "chain": "don_chain"})
                                      .merge(a_n6_single_hbond, how='inner'))
a_n6_dual_a_n3_acc_hbond_from_n6 = (a_n6_dual_a_n3_acc_res
                                    .rename(columns={"resn": "don_resn", "resi": "don_resi", "chain": "don_chain"})
                                    .merge(a_n6_dual_hbond, how='inner'))

# prepare dataframes of single and dual H-bonding A(N6) that also accept at N7 and which include H-bond information
# related to the N6 H-bond donation
a_n6_single_a_n7_acc_hbond_from_n6 = (a_n6_single_a_n7_acc_res
                                      .rename(columns={"resn": "don_resn", "resi": "don_resi", "chain": "don_chain"})
                                      .merge(a_n6_single_hbond, how='inner'))
a_n6_dual_a_n7_acc_hbond_from_n6 = (a_n6_dual_a_n7_acc_res
                                    .rename(columns={"resn": "don_resn", "resi": "don_resi", "chain": "don_chain"})
                                    .merge(a_n6_dual_hbond, how='inner'))

# prepare dataframes of single and dual H-bonding A(N6) that also accept at N1 and which include H-bond information
# related to the N1 H-bond acceptation
a_n6_single_a_n1_acc_hbond_to_n1 = (a_n6_single_a_n1_acc_res
                                    .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                    .merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]].eq(["N1", "A"])
                                           .all(axis='columns')], how='inner'))
a_n6_dual_a_n1_acc_hbond_to_n1 = (a_n6_dual_a_n1_acc_res
                                  .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                  .merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]].eq(["N1", "A"])
                                         .all(axis='columns')], how='inner'))

# prepare dataframes of single and dual H-bonding A(N6) that also accept at N3 and which include H-bond information
# related to the N3 H-bond acceptation
a_n6_single_a_n3_acc_hbond_to_n3 = (a_n6_single_a_n3_acc_res
                                    .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                    .merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]].eq(["N3", "A"])
                                           .all(axis='columns')], how='inner'))
a_n6_dual_a_n3_acc_hbond_to_n3 = (a_n6_dual_a_n3_acc_res
                                  .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                  .merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]].eq(["N3", "A"])
                                         .all(axis='columns')], how='inner'))

# prepare dataframes of single and dual H-bonding A(N6) that also accept at N7 and which include H-bond information
# related to the N7 H-bond acceptation
a_n6_single_a_n7_acc_hbond_to_n7 = (a_n6_single_a_n7_acc_res
                                    .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                    .merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]].eq(["N7", "A"])
                                           .all(axis='columns')], how='inner'))
a_n6_dual_a_n7_acc_hbond_to_n7 = (a_n6_dual_a_n7_acc_res
                                  .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                  .merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]].eq(["N7", "A"])
                                         .all(axis='columns')], how='inner'))

# prepare dataframes of single and dual H-bonding A(N6) that also accept at N1 and which include H-bond information
# related to the N6 H-bond donation and the N1 H-bond acceptation
a_n6_single_a_n1_acc_match = (a_n6_single_a_n1_acc_hbond_from_n6
                              .merge(a_n6_single_a_n1_acc_hbond_to_n1,
                                     left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
                                     right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"], how='inner'))
a_n6_dual_a_n1_acc_match = (a_n6_dual_a_n1_acc_hbond_from_n6
                            .merge(a_n6_dual_a_n1_acc_hbond_to_n1,
                                   left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
                                   right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"], how='inner'))

# prepare dataframes of single and dual H-bonding A(N6) that also accept at N3 and which include H-bond information
# related to the N6 H-bond donation and the N3 H-bond acceptation
a_n6_single_a_n3_acc_match = (a_n6_single_a_n3_acc_hbond_from_n6
                              .merge(a_n6_single_a_n3_acc_hbond_to_n3,
                                     left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
                                     right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"], how='inner'))
a_n6_dual_a_n3_acc_match = (a_n6_dual_a_n3_acc_hbond_from_n6
                            .merge(a_n6_dual_a_n3_acc_hbond_to_n3,
                                   left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
                                   right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"], how='inner'))

# prepare dataframes of single and dual H-bonding A(N6) that also accept at N7 and which include H-bond information
# related to the N6 H-bond donation and the N7 H-bond acceptation
a_n6_single_a_n7_acc_match = (a_n6_single_a_n7_acc_hbond_from_n6
                              .merge(a_n6_single_a_n7_acc_hbond_to_n7,
                                     left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
                                     right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"], how='inner'))
a_n6_dual_a_n7_acc_match = (a_n6_dual_a_n7_acc_hbond_from_n6
                            .merge(a_n6_dual_a_n7_acc_hbond_to_n7,
                                   left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
                                   right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"], how='inner'))

# Prepare dataframes with just residue information of single and dual H-bonding A(N6) that also accepts an H-bond at the
# N1. Additionally, there is overlap in the identity of the partner nucleobase, amino acid residue backbone, or amino
# acid residue side chain, where the N6 donates an H-bond and the N1 accepts an H-bond from the same partner entity.
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
                           .loc[:, ["eq_class", "don_resn_x", "don_resi_x", "don_chain_x"]]
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
                           .loc[:, ["eq_class", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())

# Prepare dataframes with just residue information of single and dual H-bonding A(N6) that also accepts an H-bond at the
# N1. Additionally, there is no overlap in the identity of the partner nucleobase, amino acid residue backbone, or amino
# acid residue side chain. In other words, the acceptor for the N6 H-bond donation and the donor for the N1 H-bond
# acceptation belong to different partner entities.
a_n6_single_a_n1_acc_no_ol = pd.concat([a_n6_single_a_n1_acc_res, a_n6_single_a_n1_acc_ol]).drop_duplicates(keep=False)
a_n6_dual_a_n1_acc_no_ol = pd.concat([a_n6_dual_a_n1_acc_res, a_n6_dual_a_n1_acc_ol]).drop_duplicates(keep=False)

print(a_n6_single_a_n1_acc_no_ol.shape[0] / a_n6_single_a_n1_acc_res.shape[0])
print(a_n6_dual_a_n1_acc_no_ol.shape[0] / a_n6_dual_a_n1_acc_res.shape[0])

# # prepare dataframes of H-bonds to A(N1) acceptors for A residues that have of no, single, and dual H-bonding A(N6)
# # a_n1_acc_a_n6_single = single_don_exo_amines_c.merge(don_hbonds_nr_c[don_hbonds_nr_c["don_resn"] == "A"],
# #                                                 on=["don_index", "eq_class"], how='inner')

# # A(N6) DONATION AND A(N1), A(N3), and A(N7) ACCEPTATION
# a_n1_acc_a_n6_no = a_n6_no_res.merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]]
#                                      .eq(["N1", "A"]).all(axis='columns')],
#                                      left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                      right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
#                                      how='inner')
# a_n1_acc_a_n6_single = a_n6_single_hbond.merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]]
#                                              .eq(["N1", "A"]).all(axis='columns')],
#                                              left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                              right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
#                                              how='inner')
# a_n1_acc_a_n6_dual = a_n6_dual_hbond.merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]]
#                                          .eq(["N1", "A"]).all(axis='columns')],
#                                          left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                          right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
#                                          how='inner')
# a_n3_acc_a_n6_no = a_n6_no_res.merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]]
#                                      .eq(["N3", "A"]).all(axis='columns')],
#                                      left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                      right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
#                                      how='inner')
# a_n3_acc_a_n6_single = a_n6_single_hbond.merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]]
#                                              .eq(["N3", "A"]).all(axis='columns')],
#                                              left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                              right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
#                                              how='inner')
# a_n3_acc_a_n6_dual = a_n6_dual_hbond.merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]]
#                                          .eq(["N3", "A"]).all(axis='columns')],
#                                          left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                          right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
#                                          how='inner')
# a_n7_acc_a_n6_no = a_n6_no_res.merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]]
#                                      .eq(["N7", "A"]).all(axis='columns')],
#                                      left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                      right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
#                                      how='inner')
# a_n7_acc_a_n6_single = a_n6_single_hbond.merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]]
#                                              .eq(["N7", "A"]).all(axis='columns')],
#                                              left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                              right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
#                                              how='inner')
# a_n7_acc_a_n6_dual = a_n6_dual_hbond.merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]]
#                                          .eq(["N7", "A"]).all(axis='columns')],
#                                          left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                          right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
#                                          how='inner')
# a_n6_donation = pd.DataFrame({
#     "Atom": ["A(N1)", "A(N1)", "A(N1)", "A(N3)", "A(N3)", "A(N3)", "A(N7)", "A(N7)", "A(N7)"],
#     "Type": ["No", "Single", "Dual", "No", "Single", "Dual", "No", "Single", "Dual"],
#     "Occurrence": [len(a_n1_acc_a_n6_no.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n1_acc_a_n6_single.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n1_acc_a_n6_dual.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n3_acc_a_n6_no.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n3_acc_a_n6_single.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n3_acc_a_n6_dual.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n7_acc_a_n6_no.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n7_acc_a_n6_single.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n7_acc_a_n6_dual.groupby(["don_index_x", "eq_class"]).groups.keys())],
#     "Total": [len(a_n6_no_res.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_single_hbond.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_dual_hbond.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_no_res.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_single_hbond.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_dual_hbond.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_no_res.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_single_hbond.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_dual_hbond.groupby(["don_index", "eq_class"]).groups.keys())]
# })
# a_n6_donation["Ratio"] = a_n6_donation["Occurrence"] / a_n6_donation["Total"]
# a_n6_donation.to_csv("plots/a_n6_donation.csv", index=False)

# prot_don_hbonds_nr = (prot_don_hbonds[prot_don_hbonds.groupby(acc_grp)["rotated_side_chain"]
#                       .transform(lambda grp: grp.str.fullmatch("none") if any(grp.str.fullmatch("none")) else True)])
                                             # .transform(lambda grp: [mem == grp.min() for mem in grp])

# A(N6) DONATION AND A(N1), A(N3), and A(N7) ACCEPTATION, DIFFERENT H-BONDING NUCLEOBASES
# a_n1_acc_a_n6_no_diff = (a_n1_acc_a_n6_no[a_n1_acc_a_n6_no.groupby(["don_index_x", "eq_class"])
#                          .transform(lambda grp: [mem["don_index_x"] == 9768 for mem in grp])])
# a_n1_acc_a_n6_no["new"] = (a_n1_acc_a_n6_no.groupby(["don_index_x", "eq_class"]).transform(lambda grp: [mem for mem in grp]))
# a_n1_acc_a_n6_dual.loc[:4, :].groupby(["don_index_x", "eq_class"]).transform(lambda grp: print(grp.loc[:, ["ang_y"]]))
# TODO set up multiple transforms to filter out groups that have one or more donor nucleobases that also accept from the exo amine

# print(a_n1_acc_a_n6_dual[(a_n1_acc_a_n6_dual["don_index_x"] == 56) & (a_n1_acc_a_n6_dual["eq_class"] == "NR_3.0_18859.1")])

# a_n1_acc_a_n6_no_diff = a_n1_acc_a_n6_no[
#     ~((a_n1_acc_a_n6_no["acc_resn_x"] == a_n1_acc_a_n6_no["don_resn_y"]) &
#       (a_n1_acc_a_n6_no["acc_resi_x"] == a_n1_acc_a_n6_no["don_resi_y"]) &
#       (a_n1_acc_a_n6_no["acc_chain_x"] == a_n1_acc_a_n6_no["don_chain_y"]) &
#       (a_n1_acc_a_n6_no["don_name_y"] != "O2'"))]
# a_n1_acc_a_n6_single_diff = a_n1_acc_a_n6_single[
#     ~((a_n1_acc_a_n6_single["acc_resn_x"] == a_n1_acc_a_n6_single["don_resn_y"]) &
#       (a_n1_acc_a_n6_single["acc_resi_x"] == a_n1_acc_a_n6_single["don_resi_y"]) &
#       (a_n1_acc_a_n6_single["acc_chain_x"] == a_n1_acc_a_n6_single["don_chain_y"]) &
#       (a_n1_acc_a_n6_single["don_name_y"] != "O2'"))]
# a_n1_acc_a_n6_dual_diff = a_n1_acc_a_n6_dual[
#     ~((a_n1_acc_a_n6_dual["acc_resn_x"] == a_n1_acc_a_n6_dual["don_resn_y"]) &
#       (a_n1_acc_a_n6_dual["acc_resi_x"] == a_n1_acc_a_n6_dual["don_resi_y"]) &
#       (a_n1_acc_a_n6_dual["acc_chain_x"] == a_n1_acc_a_n6_dual["don_chain_y"]) &
#       (a_n1_acc_a_n6_dual["don_name_y"] != "O2'"))]
# a_n3_acc_a_n6_no_diff = a_n3_acc_a_n6_no[
#     ~((a_n3_acc_a_n6_no["acc_resn_x"] == a_n3_acc_a_n6_no["don_resn_y"]) &
#       (a_n3_acc_a_n6_no["acc_resi_x"] == a_n3_acc_a_n6_no["don_resi_y"]) &
#       (a_n3_acc_a_n6_no["acc_chain_x"] == a_n3_acc_a_n6_no["don_chain_y"]) &
#       (a_n3_acc_a_n6_no["don_name_y"] != "O2'"))]
# a_n3_acc_a_n6_single_diff = a_n3_acc_a_n6_single[
#     ~((a_n3_acc_a_n6_single["acc_resn_x"] == a_n3_acc_a_n6_single["don_resn_y"]) &
#       (a_n3_acc_a_n6_single["acc_resi_x"] == a_n3_acc_a_n6_single["don_resi_y"]) &
#       (a_n3_acc_a_n6_single["acc_chain_x"] == a_n3_acc_a_n6_single["don_chain_y"]) &
#       (a_n3_acc_a_n6_single["don_name_y"] != "O2'"))]
# a_n3_acc_a_n6_dual_diff = a_n3_acc_a_n6_dual[
#     ~((a_n3_acc_a_n6_dual["acc_resn_x"] == a_n3_acc_a_n6_dual["don_resn_y"]) &
#       (a_n3_acc_a_n6_dual["acc_resi_x"] == a_n3_acc_a_n6_dual["don_resi_y"]) &
#       (a_n3_acc_a_n6_dual["acc_chain_x"] == a_n3_acc_a_n6_dual["don_chain_y"]) &
#       (a_n3_acc_a_n6_dual["don_name_y"] != "O2'"))]
# a_n7_acc_a_n6_no_diff = a_n7_acc_a_n6_no[
#     ~((a_n7_acc_a_n6_no["acc_resn_x"] == a_n7_acc_a_n6_no["don_resn_y"]) &
#       (a_n7_acc_a_n6_no["acc_resi_x"] == a_n7_acc_a_n6_no["don_resi_y"]) &
#       (a_n7_acc_a_n6_no["acc_chain_x"] == a_n7_acc_a_n6_no["don_chain_y"]) &
#       (a_n7_acc_a_n6_no["don_name_y"] != "O2'"))]
# a_n7_acc_a_n6_single_diff = a_n7_acc_a_n6_single[
#     ~((a_n7_acc_a_n6_single["acc_resn_x"] == a_n7_acc_a_n6_single["don_resn_y"]) &
#       (a_n7_acc_a_n6_single["acc_resi_x"] == a_n7_acc_a_n6_single["don_resi_y"]) &
#       (a_n7_acc_a_n6_single["acc_chain_x"] == a_n7_acc_a_n6_single["don_chain_y"]) &
#       (a_n7_acc_a_n6_single["don_name_y"] != "O2'"))]
# a_n7_acc_a_n6_dual_diff = a_n7_acc_a_n6_dual[
#     ~((a_n7_acc_a_n6_dual["acc_resn_x"] == a_n7_acc_a_n6_dual["don_resn_y"]) &
#       (a_n7_acc_a_n6_dual["acc_resi_x"] == a_n7_acc_a_n6_dual["don_resi_y"]) &
#       (a_n7_acc_a_n6_dual["acc_chain_x"] == a_n7_acc_a_n6_dual["don_chain_y"]) &
#       (a_n7_acc_a_n6_dual["don_name_y"] != "O2'"))]
# a_n6_donation_diff = pd.DataFrame({
#     "Atom": ["A(N1)", "A(N1)", "A(N1)", "A(N3)", "A(N3)", "A(N3)", "A(N7)", "A(N7)", "A(N7)"],
#     "Type": ["No", "Single", "Dual", "No", "Single", "Dual", "No", "Single", "Dual"],
#     "Occurrence": [len(a_n1_acc_a_n6_no_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n1_acc_a_n6_single_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n1_acc_a_n6_dual_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n3_acc_a_n6_no_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n3_acc_a_n6_single_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n3_acc_a_n6_dual_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n7_acc_a_n6_no_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n7_acc_a_n6_single_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n7_acc_a_n6_dual_diff.groupby(["don_index_x", "eq_class"]).groups.keys())],
#     "Total": [len(a_n6_no_res.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_single_hbond.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_dual_hbond.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_no_res.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_single_hbond.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_dual_hbond.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_no_res.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_single_hbond.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_dual_hbond.groupby(["don_index", "eq_class"]).groups.keys())]
# })
# a_n6_donation_diff["Ratio"] = a_n6_donation_diff["Occurrence"] / a_n6_donation_diff["Total"]
# a_n6_donation_diff.to_csv("plots/a_n6_donation_diff.csv", index=False)

# # A(N6) DONATION AND A(N1), A(N3), and A(N7) DONATION
# a_n1_don_a_n6_no = a_n6_no_res.merge(prot_don_hbonds_nr_c[prot_don_hbonds_nr_c[["don_name", "don_resn"]]
#                                      .eq(["N1", "A"]).all(axis='columns')],
#                                      on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                      how='inner')
# a_n1_don_a_n6_single = a_n6_single_hbond.merge(prot_don_hbonds_nr_c[prot_don_hbonds_nr_c[["don_name", "don_resn"]]
#                                              .eq(["N1", "A"]).all(axis='columns')],
#                                              on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                              how='inner')
# a_n1_don_a_n6_dual = a_n6_dual_hbond.merge(prot_don_hbonds_nr_c[prot_don_hbonds_nr_c[["don_name", "don_resn"]]
#                                          .eq(["N1", "A"]).all(axis='columns')],
#                                          on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                          how='inner')
# a_n3_don_a_n6_no = a_n6_no_res.merge(prot_don_hbonds_nr_c[prot_don_hbonds_nr_c[["don_name", "don_resn"]]
#                                      .eq(["N3", "A"]).all(axis='columns')],
#                                      on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                      how='inner')
# a_n3_don_a_n6_single = a_n6_single_hbond.merge(prot_don_hbonds_nr_c[prot_don_hbonds_nr_c[["don_name", "don_resn"]]
#                                              .eq(["N3", "A"]).all(axis='columns')],
#                                              on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                              how='inner')
# a_n3_don_a_n6_dual = a_n6_dual_hbond.merge(prot_don_hbonds_nr_c[prot_don_hbonds_nr_c[["don_name", "don_resn"]]
#                                          .eq(["N3", "A"]).all(axis='columns')],
#                                          on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                          how='inner')
# a_n7_don_a_n6_no = a_n6_no_res.merge(prot_don_hbonds_nr_c[prot_don_hbonds_nr_c[["don_name", "don_resn"]]
#                                      .eq(["N7", "A"]).all(axis='columns')],
#                                      on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                      how='inner')
# a_n7_don_a_n6_single = a_n6_single_hbond.merge(prot_don_hbonds_nr_c[prot_don_hbonds_nr_c[["don_name", "don_resn"]]
#                                              .eq(["N7", "A"]).all(axis='columns')],
#                                              on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                              how='inner')
# a_n7_don_a_n6_dual = a_n6_dual_hbond.merge(prot_don_hbonds_nr_c[prot_don_hbonds_nr_c[["don_name", "don_resn"]]
#                                          .eq(["N7", "A"]).all(axis='columns')],
#                                          on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                          how='inner')
# a_n6_donation_prot = pd.DataFrame({
#     "Atom": ["A(N1)", "A(N1)", "A(N1)", "A(N3)", "A(N3)", "A(N3)", "A(N7)", "A(N7)", "A(N7)"],
#     "Type": ["No", "Single", "Dual", "No", "Single", "Dual", "No", "Single", "Dual"],
#     "Occurrence": [len(a_n1_don_a_n6_no.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n1_don_a_n6_single.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n1_don_a_n6_dual.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n3_don_a_n6_no.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n3_don_a_n6_single.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n3_don_a_n6_dual.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n7_don_a_n6_no.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n7_don_a_n6_single.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n7_don_a_n6_dual.groupby(["don_index_x", "eq_class"]).groups.keys())],
#     "Total": [len(a_n6_no_res.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_single_hbond.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_dual_hbond.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_no_res.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_single_hbond.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_dual_hbond.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_no_res.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_single_hbond.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_dual_hbond.groupby(["don_index", "eq_class"]).groups.keys())]
# })
# a_n6_donation_prot["Ratio"] = a_n6_donation_prot["Occurrence"] / a_n6_donation_prot["Total"]
# a_n6_donation_prot.to_csv("plots/a_n6_donation_prot.csv", index=False)

# # A(N6) DONATION AND A(N1), A(N3), and A(N7) DONATION, DIFFERENT H-BONDING NUCLEOBASES
# a_n1_don_a_n6_no_diff = a_n1_don_a_n6_no[
#     ~((a_n1_don_a_n6_no["acc_resn_x"] == a_n1_don_a_n6_no["acc_resn_y"]) &
#       (a_n1_don_a_n6_no["acc_resi_x"] == a_n1_don_a_n6_no["acc_resi_y"]) &
#       (a_n1_don_a_n6_no["acc_chain_x"] == a_n1_don_a_n6_no["acc_chain_y"]) &
#       ~(a_n1_don_a_n6_no["acc_name_y"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])))]
# a_n1_don_a_n6_single_diff = a_n1_don_a_n6_single[
#     ~((a_n1_don_a_n6_single["acc_resn_x"] == a_n1_don_a_n6_single["acc_resn_y"]) &
#       (a_n1_don_a_n6_single["acc_resi_x"] == a_n1_don_a_n6_single["acc_resi_y"]) &
#       (a_n1_don_a_n6_single["acc_chain_x"] == a_n1_don_a_n6_single["acc_chain_y"]) &
#       ~(a_n1_don_a_n6_single["acc_name_y"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])))]
# a_n1_don_a_n6_dual_diff = a_n1_don_a_n6_dual[
#     ~((a_n1_don_a_n6_dual["acc_resn_x"] == a_n1_don_a_n6_dual["acc_resn_y"]) &
#       (a_n1_don_a_n6_dual["acc_resi_x"] == a_n1_don_a_n6_dual["acc_resi_y"]) &
#       (a_n1_don_a_n6_dual["acc_chain_x"] == a_n1_don_a_n6_dual["acc_chain_y"]) &
#       ~(a_n1_don_a_n6_dual["acc_name_y"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])))]
# a_n3_don_a_n6_no_diff = a_n3_don_a_n6_no[
#     ~((a_n3_don_a_n6_no["acc_resn_x"] == a_n3_don_a_n6_no["acc_resn_y"]) &
#       (a_n3_don_a_n6_no["acc_resi_x"] == a_n3_don_a_n6_no["acc_resi_y"]) &
#       (a_n3_don_a_n6_no["acc_chain_x"] == a_n3_don_a_n6_no["acc_chain_y"]) &
#       ~(a_n3_don_a_n6_no["acc_name_y"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])))]
# a_n3_don_a_n6_single_diff = a_n3_don_a_n6_single[
#     ~((a_n3_don_a_n6_single["acc_resn_x"] == a_n3_don_a_n6_single["acc_resn_y"]) &
#       (a_n3_don_a_n6_single["acc_resi_x"] == a_n3_don_a_n6_single["acc_resi_y"]) &
#       (a_n3_don_a_n6_single["acc_chain_x"] == a_n3_don_a_n6_single["acc_chain_y"]) &
#       ~(a_n3_don_a_n6_single["acc_name_y"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])))]
# a_n3_don_a_n6_dual_diff = a_n3_don_a_n6_dual[
#     ~((a_n3_don_a_n6_dual["acc_resn_x"] == a_n3_don_a_n6_dual["acc_resn_y"]) &
#       (a_n3_don_a_n6_dual["acc_resi_x"] == a_n3_don_a_n6_dual["acc_resi_y"]) &
#       (a_n3_don_a_n6_dual["acc_chain_x"] == a_n3_don_a_n6_dual["acc_chain_y"]) &
#       ~(a_n3_don_a_n6_dual["acc_name_y"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])))]
# a_n7_don_a_n6_no_diff = a_n7_don_a_n6_no[
#     ~((a_n7_don_a_n6_no["acc_resn_x"] == a_n7_don_a_n6_no["acc_resn_y"]) &
#       (a_n7_don_a_n6_no["acc_resi_x"] == a_n7_don_a_n6_no["acc_resi_y"]) &
#       (a_n7_don_a_n6_no["acc_chain_x"] == a_n7_don_a_n6_no["acc_chain_y"]) &
#       ~(a_n7_don_a_n6_no["acc_name_y"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])))]
# a_n7_don_a_n6_single_diff = a_n7_don_a_n6_single[
#     ~((a_n7_don_a_n6_single["acc_resn_x"] == a_n7_don_a_n6_single["acc_resn_y"]) &
#       (a_n7_don_a_n6_single["acc_resi_x"] == a_n7_don_a_n6_single["acc_resi_y"]) &
#       (a_n7_don_a_n6_single["acc_chain_x"] == a_n7_don_a_n6_single["acc_chain_y"]) &
#       ~(a_n7_don_a_n6_single["acc_name_y"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])))]
# a_n7_don_a_n6_dual_diff = a_n7_don_a_n6_dual[
#     ~((a_n7_don_a_n6_dual["acc_resn_x"] == a_n7_don_a_n6_dual["acc_resn_y"]) &
#       (a_n7_don_a_n6_dual["acc_resi_x"] == a_n7_don_a_n6_dual["acc_resi_y"]) &
#       (a_n7_don_a_n6_dual["acc_chain_x"] == a_n7_don_a_n6_dual["acc_chain_y"]) &
#       ~(a_n7_don_a_n6_dual["acc_name_y"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])))]
# a_n6_donation_prot_diff = pd.DataFrame({
#     "Atom": ["A(N1)", "A(N1)", "A(N1)", "A(N3)", "A(N3)", "A(N3)", "A(N7)", "A(N7)", "A(N7)"],
#     "Type": ["No", "Single", "Dual", "No", "Single", "Dual", "No", "Single", "Dual"],
#     "Occurrence": [len(a_n1_don_a_n6_no_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n1_don_a_n6_single_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n1_don_a_n6_dual_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n3_don_a_n6_no_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n3_don_a_n6_single_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n3_don_a_n6_dual_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n7_don_a_n6_no_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n7_don_a_n6_single_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                    len(a_n7_don_a_n6_dual_diff.groupby(["don_index_x", "eq_class"]).groups.keys())],
#     "Total": [len(a_n6_no_res.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_single_hbond.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_dual_hbond.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_no_res.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_single_hbond.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_dual_hbond.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_no_res.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_single_hbond.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_dual_hbond.groupby(["don_index", "eq_class"]).groups.keys())]
# })
# a_n6_donation_prot_diff["Ratio"] = a_n6_donation_prot_diff["Occurrence"] / a_n6_donation_prot_diff["Total"]
# a_n6_donation_prot_diff.to_csv("plots/a_n6_donation_prot_diff.csv", index=False)

# # H-BONDING DISTANCE, HYDROGEN VERTEX, NUCLEOBASES NOT INVOLVED IN SINGLE/DUAL DONATION
# print(a_n1_acc_a_n6_single_diff[a_n1_acc_a_n6_single_diff["vertex_y"] == "hydrogen"]["dist_y"].mean())

# TODO ignore instances when exo amine donates to sugar/phosphate of other nucleotide
# TODO all H-bond donations from exo amine must be to acceptors other than that of the nucleobase which is donating to N1
# count = a_n1_acc_a_n6_dual_diff[a_n1_acc_a_n6_dual_diff["vertex_y"] == "hydrogen"].groupby(["don_index_x", "eq_class"]).count()
# print(count[count["don_name_x"] >= 2]["don_name_x"].size)
# print(a_n1_acc_a_n6_dual_diff[(a_n1_acc_a_n6_dual_diff["don_index_x"] == 1186) & (a_n1_acc_a_n6_dual_diff["eq_class"] == "NR_3.0_41610.18")])
# print(a_n1_acc_a_n6_dual[(a_n1_acc_a_n6_dual["don_index_x"] == 56) & (a_n1_acc_a_n6_dual["eq_class"] == "NR_3.0_18859.1")])
# print(len(a_n1_acc_a_n6_dual_diff[a_n1_acc_a_n6_dual_diff["vertex_y"] == "hydrogen"].groupby(["don_index_x", "eq_class"]).groups.keys()))

# # H-BONDING DISTANCE, DONOR VERTEX, NUCLEOBASES NOT INVOLVED IN SINGLE/DUAL DONATION
# print(a_n1_acc_a_n6_single_diff[a_n1_acc_a_n6_single_diff["vertex_y"] == "donor"]["dist_y"].mean())
# print(a_n1_acc_a_n6_dual_diff[a_n1_acc_a_n6_dual_diff["vertex_y"] == "donor"]["dist_y"].mean())

# # H-BONDING DISTANCE, CANONICAL BASEPAIR
# print(au_bp_single["dist_y"].mean())
# print(au_bp_single["dist_y"].std())
# print(au_bp_dual["dist_y"].mean())
# print(au_bp_dual["dist_y"].std())

# P_VALUE
# acc = np.array([[len(a_n3_acc_a_n6_single.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                  len(a_n6_single_hbond.groupby(["don_index", "eq_class"]).groups.keys()) - len(a_n3_acc_a_n6_single.groupby(["don_index_x", "eq_class"]).groups.keys())],
#                 [len(a_n3_acc_a_n6_dual.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                  len(a_n6_dual_hbond.groupby(["don_index", "eq_class"]).groups.keys()) - len(a_n3_acc_a_n6_dual.groupby(["don_index_x", "eq_class"]).groups.keys())]])
# chi2, p_value, df, acc = scipy.stats.chi2_contingency(acc)
# print(p_value)

# TODO incorporate b-factor filter

# TODO consider revising collect_data.py to only consider nucleobases with all heavy atoms, see resn A and resi 2150 and chain A eq_class NR_3.0_90633.21 for an example of an incomplete nucleobase, leading to the following evaluation to be False
# print((a_n6_no_a_n7_acc_res["eq_class"].size + a_n6_single_a_n7_acc_res["eq_class"].size + a_n6_dual_a_n7_acc_res["eq_class"].size) == a_n7_acc_res["eq_class"].size)
