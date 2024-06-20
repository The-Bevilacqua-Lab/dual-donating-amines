"""
This script combines dataframes from individual equivalence class members, stored as csv files in the process/ folder,
into combined dataframes. The combined dataframes are then written to csv files.
"""

import pandas as pd

# Initialize combined dataframes.
h_bond_combined = pd.DataFrame(columns=["don_index", "don_name", "don_resn", "don_resi", "don_chain", "acc_index",
                                        "acc_name", "acc_resn", "acc_resi", "acc_chain", "don_acc_distance",
                                        "h_acc_distance", "h_angle", "h_dihedral", "h_name", "eq_class_members"])
nuc_combined = pd.DataFrame(columns=["don_resn", "don_resi", "don_chain", "eq_class_members"])
b_factor_combined = pd.DataFrame(columns=["resn", "resi", "chain", "subset", "mean", "eq_class_members"])
don_h_bonds_combined = pd.DataFrame(columns=["don_index", "don_name", "don_resn", "don_resi", "don_chain", "acc_index",
                                             "acc_name", "acc_resn", "acc_resi", "acc_chain", "don_acc_distance",
                                             "h_acc_distance", "h_angle", "h_dihedral", "h_name", "eq_class_members"])

# Combine the dataframes from the different equivalence class members. Do not incorporate empty dataframes.
for idx in range(len(snakemake.input.h_bond)):
    try:
        h_bond_df = pd.read_csv(snakemake.input.h_bond[idx], keep_default_na=False, na_values="NaN")
        nuc_df = pd.read_csv(snakemake.input.nuc[idx], keep_default_na=False, na_values="NaN")
        b_factor_df = pd.read_csv(snakemake.input.b_factor[idx], keep_default_na=False, na_values="NaN")
        don_h_bonds_df = pd.read_csv(snakemake.input.don_h_bonds[idx], keep_default_na=False, na_values="NaN")
        if h_bond_combined.empty and not h_bond_df.empty:
            h_bond_combined = h_bond_df.copy()
        elif not h_bond_df.empty:
            h_bond_combined = pd.concat([h_bond_combined, h_bond_df])
        if nuc_combined.empty and not nuc_df.empty:
            nuc_combined = nuc_df.copy()
        elif not nuc_df.empty:
            nuc_combined = pd.concat([nuc_combined, nuc_df])
        if b_factor_combined.empty and not b_factor_df.empty:
            b_factor_combined = b_factor_df.copy()
        elif not b_factor_df.empty:
            b_factor_combined = pd.concat([b_factor_combined, b_factor_df])
        if don_h_bonds_combined.empty and not don_h_bonds_df.empty:
            don_h_bonds_combined = don_h_bonds_df.copy()
        elif not don_h_bonds_df.empty:
            don_h_bonds_combined = pd.concat([don_h_bonds_combined, don_h_bonds_df])
    except FileNotFoundError:
        continue

# Write to csv files.
h_bond_combined.to_csv(snakemake.output.h_bond, index=False, na_rep='NaN')
nuc_combined.to_csv(snakemake.output.nuc, index=False, na_rep='NaN')
b_factor_combined.to_csv(snakemake.output.b_factor, index=False, na_rep='NaN')
don_h_bonds_combined.to_csv(snakemake.output.don_h_bonds, index=False, na_rep='NaN')
