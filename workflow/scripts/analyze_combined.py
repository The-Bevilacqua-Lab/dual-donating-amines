"""
This script will eventually do something.
"""

import os
import pandas as pd
import numpy as np
import scipy

# the working directory should house the combined/ folder
working_dir = os.getcwd()

# # create the combined dataframes
don_hbonds_nr_c = pd.read_csv(working_dir + "/combined/don_hbonds_nr_c.csv",
                              dtype={"don_resi": "object", "acc_resi": "object"})
prot_don_hbonds_nr_c = pd.read_csv(working_dir + "/combined/prot_don_hbonds_nr_c.csv",
                                   dtype={"don_resi": "object", "acc_resi": "object"})
acc_hbonds_nr_c = pd.read_csv(working_dir + "/combined/acc_hbonds_nr_c.csv",
                              dtype={"don_resi": "object", "acc_resi": "object"})
deprot_acc_hbonds_nr_c = pd.read_csv(working_dir + "/combined/deprot_acc_hbonds_nr_c.csv",
                                     dtype={"don_resi": "object", "acc_resi": "object"})
single_don_exo_amines_c = pd.read_csv(working_dir + "/combined/single_don_exo_amines_c.csv")
dual_don_exo_amines_c = pd.read_csv(working_dir + "/combined/dual_don_exo_amines_c.csv")
single_acc_carbonyls_c = pd.read_csv(working_dir + "/combined/single_acc_carbonyls_c.csv")
dual_acc_carbonyls_c = pd.read_csv(working_dir + "/combined/dual_acc_carbonyls_c.csv")
tri_acc_carbonyls_c = pd.read_csv(working_dir + "/combined/tri_acc_carbonyls_c.csv")
nuc_data_c = pd.read_csv(working_dir + "/combined/nuc_data_c.csv")

# prepare dataframes of exocyclic amines for each A, C, and G residue
a_exo_amines = nuc_data_c[nuc_data_c["atom_name"] == "N6"]
c_exo_amines = nuc_data_c[nuc_data_c["atom_name"] == "N4"]
g_exo_amines = nuc_data_c[nuc_data_c["atom_name"] == "N2"]

# prepare dataframes of carbonyls for each C and G residue
c_carbonyls = nuc_data_c[nuc_data_c["atom_name"] == "O2"]
g_carbonyls = nuc_data_c[nuc_data_c["atom_name"] == "O6"]

# prepare a df of no, single, and dual H-bonding A(N6) with all columns
a_n6_single_all = single_don_exo_amines_c.merge(don_hbonds_nr_c[don_hbonds_nr_c["don_resn"] == "A"],
                                                on=["don_index", "eq_class"], how='inner')
a_n6_dual_all = dual_don_exo_amines_c.merge(don_hbonds_nr_c[don_hbonds_nr_c["don_resn"] == "A"],
                                            on=["don_index", "eq_class"], how='inner')
a_n6_no_all = (pd.concat([a_n6_single_all,
                         a_n6_dual_all,
                         a_exo_amines.rename(columns={"index": "don_index", "atom_name": "don_name", "resn": "don_resn",
                                                      "resi": "don_resi", "chain": "don_chain"})])
               .drop_duplicates(subset=["don_index", "eq_class"], keep=False))

# prepare a list of single and dual H-bonding A residues that are involved in a canonical AU base pair
au_bp_single = (a_n6_single_all[a_n6_single_all[["acc_name", "acc_resn"]].eq(["O4", "U"]).all(axis='columns')]
                .merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["don_name", "don_resn", "acc_name"]].eq(["N3", "U", "N1"])
                       .all(axis='columns')],
                       left_on=["don_resn", "don_resi", "don_chain", "acc_resn", "acc_resi", "acc_chain", "eq_class"],
                       right_on=["acc_resn", "acc_resi", "acc_chain", "don_resn", "don_resi", "don_chain", "eq_class"],
                       how='inner'))
au_bp_dual = (a_n6_dual_all[a_n6_dual_all[["acc_name", "acc_resn"]].eq(["O4", "U"]).all(axis='columns')]
              .merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["don_name", "don_resn", "acc_name"]].eq(["N3", "U", "N1"])
                     .all(axis='columns')],
                     left_on=["don_resn", "don_resi", "don_chain", "acc_resn", "acc_resi", "acc_chain", "eq_class"],
                     right_on=["acc_resn", "acc_resi", "acc_chain", "don_resn", "don_resi", "don_chain", "eq_class"],
                     how='inner'))

# prepare a df of single and dual H-bonding C(N4) with all columns
c_n4_single_all = single_don_exo_amines_c.merge(don_hbonds_nr_c[don_hbonds_nr_c["don_resn"] == "C"],
                                                on=["don_index", "eq_class"], how='inner')
c_n4_dual_all = dual_don_exo_amines_c.merge(don_hbonds_nr_c[don_hbonds_nr_c["don_resn"] == "C"],
                                            on=["don_index", "eq_class"], how='inner')

# prepare a df of single, dual, and tri H-bonding C(O2) with all columns
c_o2_single_all = single_acc_carbonyls_c.merge(acc_hbonds_nr_c[acc_hbonds_nr_c["acc_resn"] == "C"],
                                             on=["acc_index", "eq_class"], how='inner')
c_o2_dual_all = dual_acc_carbonyls_c.merge(acc_hbonds_nr_c[acc_hbonds_nr_c["acc_resn"] == "C"],
                                             on=["acc_index", "eq_class"], how='inner')
c_o2_tri_all = tri_acc_carbonyls_c.merge(acc_hbonds_nr_c[acc_hbonds_nr_c["acc_resn"] == "C"],
                                             on=["acc_index", "eq_class"], how='inner')

# prepare a list of single and dual H-bonding C residues that are involved in a canonical CG base pair
cg_bp_single = (c_n4_single_all[c_n4_single_all[["acc_name", "acc_resn"]].eq(["O6", "G"]).all(axis='columns')]
                .merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["don_name", "don_resn", "acc_name"]].eq(["N1", "G", "N3"])
                       .all(axis='columns')],
                       left_on=["don_resn", "don_resi", "don_chain", "acc_resn", "acc_resi", "acc_chain", "eq_class"],
                       right_on=["acc_resn", "acc_resi", "acc_chain", "don_resn", "don_resi", "don_chain", "eq_class"],
                       how='inner'))
cg_bp_dual = (c_n4_dual_all[c_n4_dual_all[["acc_name", "acc_resn"]].eq(["O6", "G"]).all(axis='columns')]
              .merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["don_name", "don_resn", "acc_name"]].eq(["N1", "G", "N3"])
                     .all(axis='columns')],
                     left_on=["don_resn", "don_resi", "don_chain", "acc_resn", "acc_resi", "acc_chain", "eq_class"],
                     right_on=["acc_resn", "acc_resi", "acc_chain", "don_resn", "don_resi", "don_chain", "eq_class"],
                     how='inner'))

# prepare a df of single and dual H-bonding G(N2) with all columns
g_n2_single_all = single_don_exo_amines_c.merge(don_hbonds_nr_c[don_hbonds_nr_c["don_resn"] == "G"],
                                                on=["don_index", "eq_class"], how='inner')
g_n2_dual_all = dual_don_exo_amines_c.merge(don_hbonds_nr_c[don_hbonds_nr_c["don_resn"] == "G"],
                                            on=["don_index", "eq_class"], how='inner')

# prepare a df of single, dual, and tri H-bonding G(O6) with all columns
g_o6_single_all = single_acc_carbonyls_c.merge(acc_hbonds_nr_c[acc_hbonds_nr_c["acc_resn"] == "G"],
                                             on=["acc_index", "eq_class"], how='inner')
g_o6_dual_all = dual_acc_carbonyls_c.merge(acc_hbonds_nr_c[acc_hbonds_nr_c["acc_resn"] == "G"],
                                             on=["acc_index", "eq_class"], how='inner')
g_o6_tri_all = tri_acc_carbonyls_c.merge(acc_hbonds_nr_c[acc_hbonds_nr_c["acc_resn"] == "G"],
                                             on=["acc_index", "eq_class"], how='inner')

# prepare a list of single and dual H-bonding G residues that are involved in a canonical GC base pair
gc_bp_single = (g_n2_single_all[g_n2_single_all[["acc_name", "acc_resn"]].eq(["O2", "C"]).all(axis='columns')]
                .merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["don_name", "don_resn", "acc_name", "acc_resn"]]
                       .eq(["N1", "G", "N3", "C"]).all(axis='columns')],
                       on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
                       how='inner'))
gc_bp_dual = (g_n2_dual_all[g_n2_dual_all[["acc_name", "acc_resn"]].eq(["O2", "C"]).all(axis='columns')]
              .merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["don_name", "don_resn", "acc_name", "acc_resn"]]
                     .eq(["N1", "G", "N3", "C"]).all(axis='columns')],
                     on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
                     how='inner'))

# # write data on single and dual H-bonding and canonical base pairs to a csv file
# canonical_bp = pd.DataFrame(
#     {
#         "BP_single": [len(au_bp_single.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                       len(cg_bp_single.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                       len(gc_bp_single.groupby(["don_index_x", "eq_class"]).groups.keys())],
#         "all_single": [len(a_n6_single_all.groupby(["don_index", "eq_class"]).groups.keys()),
#                        len(c_n4_single_all.groupby(["don_index", "eq_class"]).groups.keys()),
#                        len(g_n2_single_all.groupby(["don_index", "eq_class"]).groups.keys())],
#         "BP_dual": [len(au_bp_dual.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                     len(cg_bp_dual.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                     len(gc_bp_dual.groupby(["don_index_x", "eq_class"]).groups.keys())],
#         "all_dual": [len(a_n6_dual_all.groupby(["don_index", "eq_class"]).groups.keys()),
#                      len(c_n4_dual_all.groupby(["don_index", "eq_class"]).groups.keys()),
#                      len(g_n2_dual_all.groupby(["don_index", "eq_class"]).groups.keys())]
#     },
#     index=pd.Categorical(["A", "C", "G"])
# )
# canonical_bp.to_csv("plots/canonical_bp.csv")

# # FREQUENCY OF SINGLE/MULTI EXO AMINE DONATION
# print(len(a_n6_single_all.groupby(["don_index", "eq_class"]).groups.keys()) / a_exo_amines["index"].size)
# print(len(a_n6_dual_all.groupby(["don_index", "eq_class"]).groups.keys()) / a_exo_amines["index"].size)
# print(len(c_n4_single_all.groupby(["don_index", "eq_class"]).groups.keys()) / c_exo_amines["index"].size)
# print(len(c_n4_dual_all.groupby(["don_index", "eq_class"]).groups.keys()) / c_exo_amines["index"].size)
# print(len(g_n2_single_all.groupby(["don_index", "eq_class"]).groups.keys()) / g_exo_amines["index"].size)
# print(len(g_n2_dual_all.groupby(["don_index", "eq_class"]).groups.keys()) / g_exo_amines["index"].size)

# # FREQUENCY OF SINGLE/MULTI CARBONYL ACCEPTATION
# print(len(c_o2_single_all.groupby(["acc_index", "eq_class"]).groups.keys()) / c_carbonyls["index"].size)
# print((len(c_o2_dual_all.groupby(["don_index", "eq_class"]).groups.keys()) +
#        len(c_o2_tri_all.groupby(["don_index", "eq_class"]).groups.keys())) / c_carbonyls["index"].size)
# print(len(g_o6_single_all.groupby(["acc_index", "eq_class"]).groups.keys()) / g_carbonyls["index"].size)
# print((len(g_o6_dual_all.groupby(["don_index", "eq_class"]).groups.keys()) +
#        len(g_o6_tri_all.groupby(["don_index", "eq_class"]).groups.keys())) / g_carbonyls["index"].size)

# # A(N6) DONATION AND A(N1), A(N3), and A(N7) ACCEPTATION
# a_n1_acc_a_n6_no = a_n6_no_all.merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]]
#                                      .eq(["N1", "A"]).all(axis='columns')],
#                                      left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                      right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
#                                      how='inner')
# a_n1_acc_a_n6_single = a_n6_single_all.merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]]
#                                              .eq(["N1", "A"]).all(axis='columns')],
#                                              left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                              right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
#                                              how='inner')
# a_n1_acc_a_n6_dual = a_n6_dual_all.merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]]
#                                          .eq(["N1", "A"]).all(axis='columns')],
#                                          left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                          right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
#                                          how='inner')
# a_n3_acc_a_n6_no = a_n6_no_all.merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]]
#                                      .eq(["N3", "A"]).all(axis='columns')],
#                                      left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                      right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
#                                      how='inner')
# a_n3_acc_a_n6_single = a_n6_single_all.merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]]
#                                              .eq(["N3", "A"]).all(axis='columns')],
#                                              left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                              right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
#                                              how='inner')
# a_n3_acc_a_n6_dual = a_n6_dual_all.merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]]
#                                          .eq(["N3", "A"]).all(axis='columns')],
#                                          left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                          right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
#                                          how='inner')
# a_n7_acc_a_n6_no = a_n6_no_all.merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]]
#                                      .eq(["N7", "A"]).all(axis='columns')],
#                                      left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                      right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
#                                      how='inner')
# a_n7_acc_a_n6_single = a_n6_single_all.merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]]
#                                              .eq(["N7", "A"]).all(axis='columns')],
#                                              left_on=["don_resn", "don_resi", "don_chain", "eq_class"],
#                                              right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
#                                              how='inner')
# a_n7_acc_a_n6_dual = a_n6_dual_all.merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["acc_name", "acc_resn"]]
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
#     "Total": [len(a_n6_no_all.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_single_all.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_dual_all.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_no_all.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_single_all.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_dual_all.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_no_all.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_single_all.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_dual_all.groupby(["don_index", "eq_class"]).groups.keys())]
# })
# a_n6_donation["Ratio"] = a_n6_donation["Occurrence"] / a_n6_donation["Total"]
# a_n6_donation.to_csv("plots/a_n6_donation.csv", index=False)

# # A(N6) DONATION AND A(N1), A(N3), and A(N7) ACCEPTATION, DIFFERENT H-BONDING NUCLEOBASES
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
#     "Total": [len(a_n6_no_all.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_single_all.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_dual_all.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_no_all.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_single_all.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_dual_all.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_no_all.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_single_all.groupby(["don_index", "eq_class"]).groups.keys()),
#               len(a_n6_dual_all.groupby(["don_index", "eq_class"]).groups.keys())]
# })
# a_n6_donation_diff["Ratio"] = a_n6_donation_diff["Occurrence"] / a_n6_donation_diff["Total"]
# a_n6_donation_diff.to_csv("plots/a_n6_donation_diff.csv", index=False)

# A(N6) DONATION AND A(N1), A(N3), and A(N7) DONATION
a_n1_don_a_n6_no = a_n6_no_all.merge(prot_don_hbonds_nr_c[prot_don_hbonds_nr_c[["don_name", "don_resn"]]
                                     .eq(["N1", "A"]).all(axis='columns')],
                                     on=["don_resn", "don_resi", "don_chain", "eq_class"],
                                     how='inner')
a_n1_don_a_n6_single = a_n6_single_all.merge(prot_don_hbonds_nr_c[prot_don_hbonds_nr_c[["don_name", "don_resn"]]
                                             .eq(["N1", "A"]).all(axis='columns')],
                                             on=["don_resn", "don_resi", "don_chain", "eq_class"],
                                             how='inner')
a_n1_don_a_n6_dual = a_n6_dual_all.merge(prot_don_hbonds_nr_c[prot_don_hbonds_nr_c[["don_name", "don_resn"]]
                                         .eq(["N1", "A"]).all(axis='columns')],
                                         on=["don_resn", "don_resi", "don_chain", "eq_class"],
                                         how='inner')
a_n3_don_a_n6_no = a_n6_no_all.merge(prot_don_hbonds_nr_c[prot_don_hbonds_nr_c[["don_name", "don_resn"]]
                                     .eq(["N3", "A"]).all(axis='columns')],
                                     on=["don_resn", "don_resi", "don_chain", "eq_class"],
                                     how='inner')
a_n3_don_a_n6_single = a_n6_single_all.merge(prot_don_hbonds_nr_c[prot_don_hbonds_nr_c[["don_name", "don_resn"]]
                                             .eq(["N3", "A"]).all(axis='columns')],
                                             on=["don_resn", "don_resi", "don_chain", "eq_class"],
                                             how='inner')
a_n3_don_a_n6_dual = a_n6_dual_all.merge(prot_don_hbonds_nr_c[prot_don_hbonds_nr_c[["don_name", "don_resn"]]
                                         .eq(["N3", "A"]).all(axis='columns')],
                                         on=["don_resn", "don_resi", "don_chain", "eq_class"],
                                         how='inner')
a_n7_don_a_n6_no = a_n6_no_all.merge(prot_don_hbonds_nr_c[prot_don_hbonds_nr_c[["don_name", "don_resn"]]
                                     .eq(["N7", "A"]).all(axis='columns')],
                                     on=["don_resn", "don_resi", "don_chain", "eq_class"],
                                     how='inner')
a_n7_don_a_n6_single = a_n6_single_all.merge(prot_don_hbonds_nr_c[prot_don_hbonds_nr_c[["don_name", "don_resn"]]
                                             .eq(["N7", "A"]).all(axis='columns')],
                                             on=["don_resn", "don_resi", "don_chain", "eq_class"],
                                             how='inner')
a_n7_don_a_n6_dual = a_n6_dual_all.merge(prot_don_hbonds_nr_c[prot_don_hbonds_nr_c[["don_name", "don_resn"]]
                                         .eq(["N7", "A"]).all(axis='columns')],
                                         on=["don_resn", "don_resi", "don_chain", "eq_class"],
                                         how='inner')
a_n6_donation_prot = pd.DataFrame({
    "Atom": ["A(N1)", "A(N1)", "A(N1)", "A(N3)", "A(N3)", "A(N3)", "A(N7)", "A(N7)", "A(N7)"],
    "Type": ["No", "Single", "Dual", "No", "Single", "Dual", "No", "Single", "Dual"],
    "Occurrence": [len(a_n1_don_a_n6_no.groupby(["don_index_x", "eq_class"]).groups.keys()),
                   len(a_n1_don_a_n6_single.groupby(["don_index_x", "eq_class"]).groups.keys()),
                   len(a_n1_don_a_n6_dual.groupby(["don_index_x", "eq_class"]).groups.keys()),
                   len(a_n3_don_a_n6_no.groupby(["don_index_x", "eq_class"]).groups.keys()),
                   len(a_n3_don_a_n6_single.groupby(["don_index_x", "eq_class"]).groups.keys()),
                   len(a_n3_don_a_n6_dual.groupby(["don_index_x", "eq_class"]).groups.keys()),
                   len(a_n7_don_a_n6_no.groupby(["don_index_x", "eq_class"]).groups.keys()),
                   len(a_n7_don_a_n6_single.groupby(["don_index_x", "eq_class"]).groups.keys()),
                   len(a_n7_don_a_n6_dual.groupby(["don_index_x", "eq_class"]).groups.keys())],
    "Total": [len(a_n6_no_all.groupby(["don_index", "eq_class"]).groups.keys()),
              len(a_n6_single_all.groupby(["don_index", "eq_class"]).groups.keys()),
              len(a_n6_dual_all.groupby(["don_index", "eq_class"]).groups.keys()),
              len(a_n6_no_all.groupby(["don_index", "eq_class"]).groups.keys()),
              len(a_n6_single_all.groupby(["don_index", "eq_class"]).groups.keys()),
              len(a_n6_dual_all.groupby(["don_index", "eq_class"]).groups.keys()),
              len(a_n6_no_all.groupby(["don_index", "eq_class"]).groups.keys()),
              len(a_n6_single_all.groupby(["don_index", "eq_class"]).groups.keys()),
              len(a_n6_dual_all.groupby(["don_index", "eq_class"]).groups.keys())]
})
a_n6_donation_prot["Ratio"] = a_n6_donation_prot["Occurrence"] / a_n6_donation_prot["Total"]
# a_n6_donation_prot.to_csv("plots/a_n6_donation_prot.csv", index=False)

# A(N6) DONATION AND A(N1), A(N3), and A(N7) DONATION, DIFFERENT H-BONDING NUCLEOBASES
a_n1_don_a_n6_no_diff = a_n1_don_a_n6_no[
    ~((a_n1_don_a_n6_no["acc_resn_x"] == a_n1_don_a_n6_no["acc_resn_y"]) &
      (a_n1_don_a_n6_no["acc_resi_x"] == a_n1_don_a_n6_no["acc_resi_y"]) &
      (a_n1_don_a_n6_no["acc_chain_x"] == a_n1_don_a_n6_no["acc_chain_y"]) &
      ~(a_n1_don_a_n6_no["acc_name_y"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])))]
a_n1_don_a_n6_single_diff = a_n1_don_a_n6_single[
    ~((a_n1_don_a_n6_single["acc_resn_x"] == a_n1_don_a_n6_single["acc_resn_y"]) &
      (a_n1_don_a_n6_single["acc_resi_x"] == a_n1_don_a_n6_single["acc_resi_y"]) &
      (a_n1_don_a_n6_single["acc_chain_x"] == a_n1_don_a_n6_single["acc_chain_y"]) &
      ~(a_n1_don_a_n6_single["acc_name_y"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])))]
a_n1_don_a_n6_dual_diff = a_n1_don_a_n6_dual[
    ~((a_n1_don_a_n6_dual["acc_resn_x"] == a_n1_don_a_n6_dual["acc_resn_y"]) &
      (a_n1_don_a_n6_dual["acc_resi_x"] == a_n1_don_a_n6_dual["acc_resi_y"]) &
      (a_n1_don_a_n6_dual["acc_chain_x"] == a_n1_don_a_n6_dual["acc_chain_y"]) &
      ~(a_n1_don_a_n6_dual["acc_name_y"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])))]
a_n3_don_a_n6_no_diff = a_n3_don_a_n6_no[
    ~((a_n3_don_a_n6_no["acc_resn_x"] == a_n3_don_a_n6_no["acc_resn_y"]) &
      (a_n3_don_a_n6_no["acc_resi_x"] == a_n3_don_a_n6_no["acc_resi_y"]) &
      (a_n3_don_a_n6_no["acc_chain_x"] == a_n3_don_a_n6_no["acc_chain_y"]) &
      ~(a_n3_don_a_n6_no["acc_name_y"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])))]
a_n3_don_a_n6_single_diff = a_n3_don_a_n6_single[
    ~((a_n3_don_a_n6_single["acc_resn_x"] == a_n3_don_a_n6_single["acc_resn_y"]) &
      (a_n3_don_a_n6_single["acc_resi_x"] == a_n3_don_a_n6_single["acc_resi_y"]) &
      (a_n3_don_a_n6_single["acc_chain_x"] == a_n3_don_a_n6_single["acc_chain_y"]) &
      ~(a_n3_don_a_n6_single["acc_name_y"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])))]
a_n3_don_a_n6_dual_diff = a_n3_don_a_n6_dual[
    ~((a_n3_don_a_n6_dual["acc_resn_x"] == a_n3_don_a_n6_dual["acc_resn_y"]) &
      (a_n3_don_a_n6_dual["acc_resi_x"] == a_n3_don_a_n6_dual["acc_resi_y"]) &
      (a_n3_don_a_n6_dual["acc_chain_x"] == a_n3_don_a_n6_dual["acc_chain_y"]) &
      ~(a_n3_don_a_n6_dual["acc_name_y"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])))]
a_n7_don_a_n6_no_diff = a_n7_don_a_n6_no[
    ~((a_n7_don_a_n6_no["acc_resn_x"] == a_n7_don_a_n6_no["acc_resn_y"]) &
      (a_n7_don_a_n6_no["acc_resi_x"] == a_n7_don_a_n6_no["acc_resi_y"]) &
      (a_n7_don_a_n6_no["acc_chain_x"] == a_n7_don_a_n6_no["acc_chain_y"]) &
      ~(a_n7_don_a_n6_no["acc_name_y"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])))]
a_n7_don_a_n6_single_diff = a_n7_don_a_n6_single[
    ~((a_n7_don_a_n6_single["acc_resn_x"] == a_n7_don_a_n6_single["acc_resn_y"]) &
      (a_n7_don_a_n6_single["acc_resi_x"] == a_n7_don_a_n6_single["acc_resi_y"]) &
      (a_n7_don_a_n6_single["acc_chain_x"] == a_n7_don_a_n6_single["acc_chain_y"]) &
      ~(a_n7_don_a_n6_single["acc_name_y"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])))]
a_n7_don_a_n6_dual_diff = a_n7_don_a_n6_dual[
    ~((a_n7_don_a_n6_dual["acc_resn_x"] == a_n7_don_a_n6_dual["acc_resn_y"]) &
      (a_n7_don_a_n6_dual["acc_resi_x"] == a_n7_don_a_n6_dual["acc_resi_y"]) &
      (a_n7_don_a_n6_dual["acc_chain_x"] == a_n7_don_a_n6_dual["acc_chain_y"]) &
      ~(a_n7_don_a_n6_dual["acc_name_y"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])))]
a_n6_donation_prot_diff = pd.DataFrame({
    "Atom": ["A(N1)", "A(N1)", "A(N1)", "A(N3)", "A(N3)", "A(N3)", "A(N7)", "A(N7)", "A(N7)"],
    "Type": ["No", "Single", "Dual", "No", "Single", "Dual", "No", "Single", "Dual"],
    "Occurrence": [len(a_n1_don_a_n6_no_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
                   len(a_n1_don_a_n6_single_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
                   len(a_n1_don_a_n6_dual_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
                   len(a_n3_don_a_n6_no_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
                   len(a_n3_don_a_n6_single_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
                   len(a_n3_don_a_n6_dual_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
                   len(a_n7_don_a_n6_no_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
                   len(a_n7_don_a_n6_single_diff.groupby(["don_index_x", "eq_class"]).groups.keys()),
                   len(a_n7_don_a_n6_dual_diff.groupby(["don_index_x", "eq_class"]).groups.keys())],
    "Total": [len(a_n6_no_all.groupby(["don_index", "eq_class"]).groups.keys()),
              len(a_n6_single_all.groupby(["don_index", "eq_class"]).groups.keys()),
              len(a_n6_dual_all.groupby(["don_index", "eq_class"]).groups.keys()),
              len(a_n6_no_all.groupby(["don_index", "eq_class"]).groups.keys()),
              len(a_n6_single_all.groupby(["don_index", "eq_class"]).groups.keys()),
              len(a_n6_dual_all.groupby(["don_index", "eq_class"]).groups.keys()),
              len(a_n6_no_all.groupby(["don_index", "eq_class"]).groups.keys()),
              len(a_n6_single_all.groupby(["don_index", "eq_class"]).groups.keys()),
              len(a_n6_dual_all.groupby(["don_index", "eq_class"]).groups.keys())]
})
a_n6_donation_prot_diff["Ratio"] = a_n6_donation_prot_diff["Occurrence"] / a_n6_donation_prot_diff["Total"]
# a_n6_donation_prot_diff.to_csv("plots/a_n6_donation_prot_diff.csv", index=False)

# # H-BONDING DISTANCE, NUCLEOBASES NOT INVOLVED IN SINGLE/DUAL DONATION
# print(a_n1_acc_a_n6_single_2[a_n1_acc_a_n6_single_2["vertex_y"] == "hydrogen"]["dist_y"].mean())
# print(a_n1_acc_a_n6_dual_2[a_n1_acc_a_n6_dual_2["vertex_y"] == "hydrogen"]["dist_y"].mean())
# print(a_n1_acc_a_n6_single_2[a_n1_acc_a_n6_single_2["vertex_y"] == "donor"]["dist_y"].mean())
# print(a_n1_acc_a_n6_dual_2[a_n1_acc_a_n6_dual_2["vertex_y"] == "donor"]["dist_y"].mean())

# # H-BONDING DISTANCE, CANONICAL BASEPAIR
# print(au_bp_single["dist_y"].mean())
# print(au_bp_single["dist_y"].std())
# print(au_bp_dual["dist_y"].mean())
# print(au_bp_dual["dist_y"].std())
# print(cg_bp_single["dist_y"].mean())
# print(cg_bp_single["dist_y"].std())
# print(cg_bp_dual["dist_y"].mean())
# print(cg_bp_dual["dist_y"].std())
# print(gc_bp_single["dist_y"].mean())
# print(gc_bp_single["dist_y"].std())
# print(gc_bp_dual["dist_y"].mean())
# print(gc_bp_dual["dist_y"].std())

# P_VALUE
# acc = np.array([[len(a_n3_acc_a_n6_single.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                  len(a_n6_single_all.groupby(["don_index", "eq_class"]).groups.keys()) - len(a_n3_acc_a_n6_single.groupby(["don_index_x", "eq_class"]).groups.keys())],
#                 [len(a_n3_acc_a_n6_dual.groupby(["don_index_x", "eq_class"]).groups.keys()),
#                  len(a_n6_dual_all.groupby(["don_index", "eq_class"]).groups.keys()) - len(a_n3_acc_a_n6_dual.groupby(["don_index_x", "eq_class"]).groups.keys())]])
# chi2, p_value, df, acc = scipy.stats.chi2_contingency(acc)
# print(p_value)

# TODO incorporate b-factor filter
