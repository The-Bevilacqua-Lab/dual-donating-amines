"""
This script organizes the data and applies the H-bonding criteria.
"""

import os
import pandas as pd

# Set the name of the H-bond file.
h_bond_file = snakemake.input.h_bond

# Extract the data from the count csv file.
if os.path.getsize(snakemake.input.counts) > 0:
    count_data = pd.read_csv(snakemake.input.counts, comment="#", na_filter=False, dtype={"chain": "object"})
else:
    count_data = pd.DataFrame(columns=["index", "name", "resn", "resi", "chain", "count_1", "count_2"])

# Extract the data from the h_bond csv file and remove redundant lines.
if os.path.getsize(h_bond_file) > 0:
    unique_col_comb = ["don_index", "acc_index", "h_name"]
    # noinspection PyTypeChecker
    h_bond_data = (pd.read_csv(h_bond_file, comment="#", keep_default_na=False,
                               na_values={"h_acc_distance": "NaN", "h_angle": "NaN", "h_dihedral": "NaN"},
                               dtype={"don_chain": "object", "acc_chain": "object"})
                   .drop_duplicates(subset=unique_col_comb))
else:
    h_bond_data = pd.DataFrame(columns=["don_index", "don_name", "don_resn", "don_resi", "don_chain", "acc_index",
                                        "acc_name", "acc_resn", "acc_resi", "acc_chain", "don_acc_distance",
                                        "h_acc_distance", "h_angle", "h_dihedral", "h_name"])

# Extract the data from the nuc csv file.
if os.path.getsize(snakemake.input.nuc) > 0:
    nuc_data = pd.read_csv(snakemake.input.nuc, comment="#", na_filter=False, dtype={"chain": "object"})
else:
    nuc_data = pd.DataFrame(columns=["don_resn", "don_resi", "don_chain"])

# Create a dictionary based on the file containing information on nucleobase b-factors.
if os.path.getsize(snakemake.input.b_factor) > 0:
    b_factor_grp = ["resn", "resi", "chain", "subset"]
    b_factor_data = pd.read_csv(snakemake.input.b_factor, comment="#", na_filter=False, dtype={"chain": "object"})
    b_factor_data["mean"] = b_factor_data.groupby(b_factor_grp)["b-factor"].transform("mean")
    b_factor_data = b_factor_data.drop_duplicates(subset=b_factor_grp).drop(columns=["index", "name", "b-factor"])
else:
    b_factor_data = pd.DataFrame(columns=["resn", "resi", "chain", "subset", "mean"])

# Specify the ID of the equivalence class member within an additional column in the dataframes.
count_data["eq_class_members"] = snakemake.wildcards.eq_class_members
h_bond_data["eq_class_members"] = snakemake.wildcards.eq_class_members
nuc_data["eq_class_members"] = snakemake.wildcards.eq_class_members
b_factor_data["eq_class_members"] = snakemake.wildcards.eq_class_members

# Write to csv files.
count_data.to_csv(snakemake.output.counts, index=False, na_rep='NaN')
h_bond_data.to_csv(snakemake.output.h_bond, index=False, na_rep='NaN')
nuc_data.to_csv(snakemake.output.nuc, index=False, na_rep='NaN')
b_factor_data.to_csv(snakemake.output.b_factor, index=False, na_rep='NaN')
