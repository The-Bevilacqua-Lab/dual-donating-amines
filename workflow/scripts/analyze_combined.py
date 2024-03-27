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
prot_don_hbonds_c = pd.read_csv(working_dir + "/combined/prot_don_hbonds_c.csv")
acc_hbonds_nr_c = pd.read_csv(working_dir + "/combined/acc_hbonds_nr_c.csv")
deprot_acc_hbonds_c = pd.read_csv(working_dir + "/combined/deprot_acc_hbonds_c.csv")
single_don_exo_amines_c = pd.read_csv(working_dir + "/combined/single_don_exo_amines_c.csv")
dual_don_exo_amines_c = pd.read_csv(working_dir + "/combined/dual_don_exo_amines_c.csv")
# don_hbonds_geom_c = pd.read_csv(working_dir + "/combined/don_hbonds_geom_c.csv")
# don_hbonds_geom_filtered_c = pd.read_csv(working_dir + "/combined/don_hbonds_geom_filtered_c.csv")
nuc_data_c = pd.read_csv(working_dir + "/combined/nuc_data_c.csv")

exo_amines = nuc_data_c[nuc_data_c["atom_name"].isin(["N2", "N4", "N6"])]

# TODO the acceptor index and hydrogen needs to be specified in the single and dual donor exo amine variables in process_data.py
single_all = single_don_exo_amines_c.merge(don_hbonds_nr_c, on=["don_index", "eq_class"], how='inner')
dual_all = dual_don_exo_amines_c.merge(don_hbonds_nr_c, on=["don_index", "eq_class"], how='inner')

# # PROT DON
# prot_don_single = single_all.merge(prot_don_hbonds_c, on=["don_resn", "don_resi", "don_chain", "eq_class"], how='inner')
# prot_don_dual = dual_all.merge(prot_don_hbonds_c, on=["don_resn", "don_resi", "don_chain", "eq_class"], how='inner')
# # print(prot_don_single["don_resn"].size/single_don_exo_amines_c["don_index"].size)
# # print(prot_don_dual["don_resn"].size/dual_don_exo_amines_c["don_index"].size)
# prot_don = np.array([[prot_don_single["don_resn"].size, single_don_exo_amines_c["don_index"].size - prot_don_single["don_resn"].size], [prot_don_dual["don_resn"].size, dual_don_exo_amines_c["don_index"].size - prot_don_dual["don_resn"].size]])
# chi2, p_value, df, prot_don = scipy.stats.chi2_contingency(prot_don)
# print(p_value)

# # ACC
# acc_single = single_all.merge(acc_hbonds_nr_c, left_on=["don_resn", "don_resi", "don_chain", "eq_class"], right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"], how='inner')
# acc_dual = dual_all.merge(acc_hbonds_nr_c, left_on=["don_resn", "don_resi", "don_chain", "eq_class"], right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"], how='inner')
# # print(acc_single["don_resn_x"].size/single_don_exo_amines_c["don_index"].size)
# # print(acc_dual["don_resn_x"].size/dual_don_exo_amines_c["don_index"].size)
# acc = np.array([[acc_single["don_resn_x"].size, single_don_exo_amines_c["don_index"].size*8 - acc_single["don_resn_x"].size], [acc_dual["don_resn_x"].size, dual_don_exo_amines_c["don_index"].size*8 - acc_dual["don_resn_x"].size]])
# chi2, p_value, df, acc = scipy.stats.chi2_contingency(acc)
# print(p_value)

# DEPROT_ACC
deprot_acc_single = single_all.merge(deprot_acc_hbonds_c, left_on=["don_resn", "don_resi", "don_chain", "eq_class"], right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"], how='inner')
deprot_acc_dual = dual_all.merge(deprot_acc_hbonds_c, left_on=["don_resn", "don_resi", "don_chain", "eq_class"], right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class"], how='inner')
# print(deprot_acc_single["don_resn_x"].size/single_don_exo_amines_c["don_index"].size)
# print(deprot_acc_dual["don_resn_x"].size/dual_don_exo_amines_c["don_index"].size)
deprot_acc = np.array([[deprot_acc_single["don_resn_x"].size, single_all[single_all["don_resn"] == "G"]["don_resn"].size - deprot_acc_single["don_resn_x"].size], [deprot_acc_dual["don_resn_x"].size, single_all[single_all["don_resn"] == "G"]["don_resn"].size - deprot_acc_dual["don_resn_x"].size]])
chi2, p_value, df, deprot_acc = scipy.stats.chi2_contingency(deprot_acc)
print(p_value)

# unique_col_comb = ["index", "eq_class"]
# new = df.drop_duplicates(subset=unique_col_comb)
# .merge(pd.DataFrame(const.DONORS_OF_INTEREST, columns=["don_resn", "don_name"]), how='inner'))
