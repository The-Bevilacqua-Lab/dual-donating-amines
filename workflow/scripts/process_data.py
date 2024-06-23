"""
This script organizes the data and applies the H-bonding criteria.
"""

import pandas as pd
import residue_library

# Set the H-bonding criteria.
H_DIST_MAX = snakemake.config["h_dist_max"]
H_ANG_MIN = snakemake.config["h_ang_min"]

# Set the name of the H-bond file.
h_bond_file = snakemake.input.h_bond

# Extract the data from the count csv file.
count_data = pd.read_csv(snakemake.input.count, comment="#", na_filter=False, dtype={"chain": "object"})

# Extract the data from the h_bond csv file and remove redundant lines.
unique_col_comb = ["don_index", "acc_index", "h_name"]
# noinspection PyTypeChecker
h_bond_data = (pd.read_csv(h_bond_file, comment="#", keep_default_na=False,
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

# Prepare a list of chains based on the string describing the equivalence class member.
chain_build = []
index = 0
track = 0
for char in snakemake.wildcards.eq_class_members:
    if char != '_' and index > 4:
        if index != track:
            chain_build.append([char])
            track = index
        else:
            chain_build[-1].append(char)
    if char == '_':
        index += 1
chain_list = []
for chain in chain_build:
    chain_list.append("".join(chain))

# Prepare a list containing the residue and atom names of the donors of interest.
donors_of_interest = []
for donor in snakemake.config["donors_of_interest"]:
    donors_of_interest.append([donor.split(".")[0], donor.split(".")[1]])

# Identify atom pairs that meet the H-bond criteria and include a donor of interest.
don_h_bonds = (h_bond_data[(h_bond_data["h_acc_distance"] <= H_DIST_MAX) &
                           (h_bond_data["h_angle"] >= H_ANG_MIN)]
               .merge(pd.DataFrame(donors_of_interest, columns=["don_resn", "don_name"]), how='inner')
               .merge(pd.DataFrame({'don_chain': chain_list}), how='inner'))

# Specify the ID of the equivalence class member within an additional column in the dataframes.
count_data["eq_class_members"] = snakemake.wildcards.eq_class_members
h_bond_data["eq_class_members"] = snakemake.wildcards.eq_class_members
nuc_data["eq_class_members"] = snakemake.wildcards.eq_class_members
b_factor_data["eq_class_members"] = snakemake.wildcards.eq_class_members
don_h_bonds["eq_class_members"] = snakemake.wildcards.eq_class_members

# Write to csv files.
count_data.to_csv(snakemake.output.count, index=False, na_rep='NaN')
h_bond_data.to_csv(snakemake.output.h_bond, index=False, na_rep='NaN')
nuc_data.to_csv(snakemake.output.nuc, index=False, na_rep='NaN')
b_factor_data.to_csv(snakemake.output.b_factor, index=False, na_rep='NaN')
don_h_bonds.to_csv(snakemake.output.don_h_bonds, index=False, na_rep='NaN')
