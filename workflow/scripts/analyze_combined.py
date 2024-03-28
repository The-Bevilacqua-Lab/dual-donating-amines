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
don_hbonds_nr_c = pd.read_csv(working_dir + "/combined/don_hbonds_nr_c.csv")
# prot_don_hbonds_c = pd.read_csv(working_dir + "/combined/prot_don_hbonds_c.csv")
acc_hbonds_nr_c = pd.read_csv(working_dir + "/combined/acc_hbonds_nr_c.csv")
# deprot_acc_hbonds_c = pd.read_csv(working_dir + "/combined/deprot_acc_hbonds_c.csv")
single_don_exo_amines_c = pd.read_csv(working_dir + "/combined/single_don_exo_amines_c.csv")
dual_don_exo_amines_c = pd.read_csv(working_dir + "/combined/dual_don_exo_amines_c.csv")
# don_hbonds_geom_c = pd.read_csv(working_dir + "/combined/don_hbonds_geom_c.csv")
# don_hbonds_geom_filtered_c = pd.read_csv(working_dir + "/combined/don_hbonds_geom_filtered_c.csv")
# nuc_data_c = pd.read_csv(working_dir + "/combined/nuc_data_c.csv")

# # prepare a df of all residues that contain an exocyclic amine
# exo_amines = nuc_data_c[nuc_data_c["atom_name"].isin(["N2", "N4", "N6"])]

# prepare a df of single and dual H-bonding A residues with all columns
a_single_all = single_don_exo_amines_c.merge(don_hbonds_nr_c[don_hbonds_nr_c["don_resn"] == "A"],
                                             on=["don_index", "eq_class"], how='inner')
a_dual_all = dual_don_exo_amines_c.merge(don_hbonds_nr_c[don_hbonds_nr_c["don_resn"] == "A"],
                                         on=["don_index", "eq_class"], how='inner')

# prepare a list of single and dual H-bonding A residues that are involved in a canonical AU base pair
au_bp_single = (a_single_all[a_single_all[["acc_name", "acc_resn"]].eq(["O4", "U"]).all(axis='columns')]
                .merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["don_name", "don_resn", "acc_name"]].eq(["N3", "U", "N1"])
                       .all(axis='columns')],
                       left_on=["don_resn", "don_resi", "don_chain", "acc_resn", "acc_resi", "acc_chain", "eq_class"],
                       right_on=["acc_resn", "acc_resi", "acc_chain", "don_resn", "don_resi", "don_chain", "eq_class"],
                       how='inner'))
au_bp_dual = (a_dual_all[a_dual_all[["acc_name", "acc_resn"]].eq(["O4", "U"]).all(axis='columns')]
              .merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["don_name", "don_resn", "acc_name"]].eq(["N3", "U", "N1"])
                     .all(axis='columns')],
                     left_on=["don_resn", "don_resi", "don_chain", "acc_resn", "acc_resi", "acc_chain", "eq_class"],
                     right_on=["acc_resn", "acc_resi", "acc_chain", "don_resn", "don_resi", "don_chain", "eq_class"],
                     how='inner'))

# prepare a df of single and dual H-bonding C residues with all columns
c_single_all = single_don_exo_amines_c.merge(don_hbonds_nr_c[don_hbonds_nr_c["don_resn"] == "C"],
                                             on=["don_index", "eq_class"], how='inner')
c_dual_all = dual_don_exo_amines_c.merge(don_hbonds_nr_c[don_hbonds_nr_c["don_resn"] == "C"],
                                         on=["don_index", "eq_class"], how='inner')

# prepare a list of single and dual H-bonding C residues that are involved in a canonical CG base pair
cg_bp_single = (c_single_all[c_single_all[["acc_name", "acc_resn"]].eq(["O6", "G"]).all(axis='columns')]
                .merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["don_name", "don_resn", "acc_name"]].eq(["N1", "G", "N3"])
                       .all(axis='columns')],
                       left_on=["don_resn", "don_resi", "don_chain", "acc_resn", "acc_resi", "acc_chain", "eq_class"],
                       right_on=["acc_resn", "acc_resi", "acc_chain", "don_resn", "don_resi", "don_chain", "eq_class"],
                       how='inner'))
cg_bp_dual = (c_dual_all[c_dual_all[["acc_name", "acc_resn"]].eq(["O6", "G"]).all(axis='columns')]
              .merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["don_name", "don_resn", "acc_name"]].eq(["N1", "G", "N3"])
                     .all(axis='columns')],
                     left_on=["don_resn", "don_resi", "don_chain", "acc_resn", "acc_resi", "acc_chain", "eq_class"],
                     right_on=["acc_resn", "acc_resi", "acc_chain", "don_resn", "don_resi", "don_chain", "eq_class"],
                     how='inner'))

# prepare a df of single and dual H-bonding G residues with all columns
g_single_all = single_don_exo_amines_c.merge(don_hbonds_nr_c[don_hbonds_nr_c["don_resn"] == "G"],
                                             on=["don_index", "eq_class"], how='inner')
g_dual_all = dual_don_exo_amines_c.merge(don_hbonds_nr_c[don_hbonds_nr_c["don_resn"] == "G"],
                                         on=["don_index", "eq_class"], how='inner')

# prepare a list of single and dual H-bonding G residues that are involved in a canonical GC base pair
gc_bp_single = (g_single_all[g_single_all[["acc_name", "acc_resn"]].eq(["O2", "C"]).all(axis='columns')]
                .merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["don_name", "don_resn", "acc_name", "acc_resn"]]
                       .eq(["N1", "G", "N3", "C"]).all(axis='columns')],
                       on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
                       how='inner'))
gc_bp_dual = (g_dual_all[g_dual_all[["acc_name", "acc_resn"]].eq(["O2", "C"]).all(axis='columns')]
              .merge(acc_hbonds_nr_c[acc_hbonds_nr_c[["don_name", "don_resn", "acc_name", "acc_resn"]]
                     .eq(["N1", "G", "N3", "C"]).all(axis='columns')],
                     on=["acc_resn", "acc_resi", "acc_chain", "eq_class"],
                     how='inner'))

# # write data on single and dual H-bonding and canonical base pairs to a csv file
# canonical_bp = pd.DataFrame(
#     {
#         "BP_single": [au_bp_single["eq_class"].size, cg_bp_single["eq_class"].size, gc_bp_single["eq_class"].size],
#         "all_single": [a_single_all["eq_class"].size, c_single_all["eq_class"].size, g_single_all["eq_class"].size],
#         "BP_dual": [au_bp_dual["eq_class"].size, cg_bp_dual["eq_class"].size, gc_bp_dual["eq_class"].size],
#         "all_dual": [len(a_dual_all.groupby(["don_index", "eq_class"]).groups.keys()), len(c_dual_all.groupby(["don_index", "eq_class"]).groups.keys()), len(g_dual_all.groupby(["don_index", "eq_class"]).groups.keys())]
#     },
#     index=pd.Categorical(["A", "C", "G"])
# )
# canonical_bp.to_csv("plots/canonical_bp.csv")

# # P_VALUE PROT DON
# # TODO take into account that every two entries in dual_all accounts for one residue
# prot_don_single = single_all.merge(prot_don_hbonds_c, on=["don_resn", "don_resi", "don_chain", "eq_class"], how='inner')
# prot_don_dual = dual_all.merge(prot_don_hbonds_c, on=["don_resn", "don_resi", "don_chain", "eq_class"], how='inner')
# # print(prot_don_single["don_resn"].size/single_don_exo_amines_c["don_index"].size)
# # print(prot_don_dual["don_resn"].size/dual_don_exo_amines_c["don_index"].size)
# prot_don = np.array([[prot_don_single["don_resn"].size, single_don_exo_amines_c["don_index"].size - prot_don_single["don_resn"].size], [prot_don_dual["don_resn"].size, dual_don_exo_amines_c["don_index"].size - prot_don_dual["don_resn"].size]])
# chi2, p_value, df, prot_don = scipy.stats.chi2_contingency(prot_don)
# print(p_value)

# # P_VALUE ACC
# # TODO take into account that every two entries in dual_all accounts for one residue
# acc_single = single_all.merge(acc_hbonds_nr_c, left_on=["don_resn", "don_resi", "don_chain", "eq_class"], right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"], how='inner')
# acc_dual = dual_all.merge(acc_hbonds_nr_c, left_on=["don_resn", "don_resi", "don_chain", "eq_class"], right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"], how='inner')
# # print(acc_single["don_resn_x"].size/single_don_exo_amines_c["don_index"].size)
# # print(acc_dual["don_resn_x"].size/dual_don_exo_amines_c["don_index"].size)
# acc = np.array([[acc_single["don_resn_x"].size, single_don_exo_amines_c["don_index"].size*8 - acc_single["don_resn_x"].size], [acc_dual["don_resn_x"].size, dual_don_exo_amines_c["don_index"].size*8 - acc_dual["don_resn_x"].size]])
# chi2, p_value, df, acc = scipy.stats.chi2_contingency(acc)
# print(p_value)

# # P_VALUE DEPROT_ACC
# # TODO take into account that every two entries in dual_all accounts for one residue
# deprot_acc_single = single_all.merge(deprot_acc_hbonds_c, left_on=["don_resn", "don_resi", "don_chain", "eq_class"], right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"], how='inner')
# deprot_acc_dual = dual_all.merge(deprot_acc_hbonds_c, left_on=["don_resn", "don_resi", "don_chain", "eq_class"], right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"], how='inner')
# # print(deprot_acc_single["don_resn_x"].size/single_don_exo_amines_c["don_index"].size)
# # print(deprot_acc_dual["don_resn_x"].size/dual_don_exo_amines_c["don_index"].size)
# deprot_acc = np.array([[deprot_acc_single["don_resn_x"].size, single_all[single_all["don_resn"] == "G"]["don_resn"].size - deprot_acc_single["don_resn_x"].size], [deprot_acc_dual["don_resn_x"].size, single_all[single_all["don_resn"] == "G"]["don_resn"].size - deprot_acc_dual["don_resn_x"].size]])
# chi2, p_value, df, deprot_acc = scipy.stats.chi2_contingency(deprot_acc)
# print(p_value)

# unique_col_comb = ["index", "eq_class"]
# new = df.drop_duplicates(subset=unique_col_comb)
# .merge(pd.DataFrame(const.DONORS_OF_INTEREST, columns=["don_resn", "don_name"]), how='inner'))

# TODO incorporate b-factor filter
