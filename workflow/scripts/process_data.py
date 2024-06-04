"""
This module will eventually do something.
"""

import os
import pandas as pd
import const

# The working directory should house the data/ folder.
working_dir = os.getcwd()

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

# Extract the data from the equivalence class file.
eq_class_members_data = (pd.read_csv(snakemake.input.eq_class_members, header=None, comment="#",
                                     skiprows=[3 if snakemake.config["commit_hash"] else 2],
                                     na_filter=False, dtype="object")
                         .set_index(pd.Series(["pdb_id", "model", "chain"])).transpose())

# Identify atom pairs that meet the H-bond criteria and include a donor of interest.
don_h_bonds = (h_bond_data[(h_bond_data["h_acc_distance"] <= H_DIST_MAX) &
                           (h_bond_data["h_angle"] >= 180.0 - H_ANG_TOL)]
               .merge(pd.DataFrame(const.DONORS_OF_INTEREST, columns=["don_resn", "don_name"]), how='inner')
               .merge(eq_class_members_data, left_on="don_chain", right_on="chain", how='inner')
               .drop(columns=["pdb_id", "model", "chain"]))

# Identify atom pairs that meet the H-bond criteria and include a protonated donor of interest. Values in the _merge
# column matching left_only indicate acceptors that cannot typically also donate an H-bond.
prot_don_h_bonds_unfiltered = (h_bond_data[(h_bond_data["h_acc_distance"] <= H_DIST_MAX) &
                                           (h_bond_data["h_angle"] >= 180.0 - H_ANG_TOL)]
                               .merge(pd.DataFrame(const.PROT_DONORS_OF_INTEREST, columns=["don_resn", "don_name"]),
                                      how='inner')
                               .merge(pd.DataFrame(don_acc_atoms, columns=["acc_resn", "acc_name"]), how='left',
                                      indicator=True))

# Filter prot_don_h_bonds_unfiltered such that only acceptors that cannot typically also donate an H-bond are included.
prot_don_h_bonds = (prot_don_h_bonds_unfiltered[prot_don_h_bonds_unfiltered["_merge"] == "left_only"]
                    .drop(columns="_merge")
                    .merge(eq_class_members_data, left_on="don_chain", right_on="chain", how='inner')
                    .drop(columns=["pdb_id", "model", "chain"]))

# Identify atom pairs that meet the H-bond criteria and include an acceptor of interest.
acc_h_bonds = (h_bond_data[(h_bond_data["h_acc_distance"].notna()) & (h_bond_data["h_acc_distance"] <= H_DIST_MAX) &
                           (h_bond_data["h_angle"] >= 180.0 - H_ANG_TOL)]
               .merge(pd.DataFrame(const.ACCEPTORS_OF_INTEREST, columns=["acc_resn", "acc_name"]), how='inner')
               .merge(eq_class_members_data, left_on="acc_chain", right_on="chain", how='inner')
               .drop(columns=["pdb_id", "model", "chain"]))

# Write dataframes to csv files.
nuc_data.assign(eq_class_members=snakemake.wildcards.eq_class_members).to_csv(snakemake.output.nuc_data, index=False)
don_h_bonds.assign(eq_class_members=snakemake.wildcards.eq_class_members).to_csv(snakemake.output.don_h_bonds_nr, index=False)
prot_don_h_bonds.assign(eq_class_members=snakemake.wildcards.eq_class_members).to_csv(snakemake.output.prot_don_h_bonds_nr, index=False)
acc_h_bonds.assign(eq_class_members=snakemake.wildcards.eq_class_members).to_csv(snakemake.output.acc_h_bonds_nr, index=False)
b_factor_data.assign(eq_class_members=snakemake.wildcards.eq_class_members).to_csv(snakemake.output.b_factor_data, index=False)
h_bond_data.assign(eq_class_members=snakemake.wildcards.eq_class_members).to_csv(snakemake.output.h_bond_data, index=False)

# Create a dictionary that will contain all the data categories and their associated files.
categories = {
    "nuc_data": [],
    "b_factor_data": [],
    "don_hbonds_nr": [],
    "prot_don_hbonds_nr": [],
    "acc_hbonds_nr": [],
    "hbond_data": []
}

# Prepare a list of files from each category folder to combine.
for cat in categories.keys():
    for file in os.listdir(working_dir + "/analysis/" + cat):
        # select files with the correct suffix
        if file[-(len(cat) + 5):] == "_" + cat + ".csv":
            # find out which files actually contain data
            contains_data = False
            with open(working_dir + "/analysis/" + cat + "/" + file, "r") as read_file:
                num_lines = 0
                for lines in read_file:
                    if num_lines > 0:
                        contains_data = True
                        break
                    else:
                        num_lines += 1
            # only include files that contain data
            if contains_data:
                categories[cat].append(working_dir + "/analysis/" + cat + "/" + file)

# Create the combined/ folder to write csv files to.
os.mkdir(working_dir + "/combined")

# Create the combined dataframes and write to csv files.
pd.concat([pd.read_csv(file, na_filter=False) for file in categories["nuc_data"]]).to_csv("combined/nuc_data.csv", index=False)
pd.concat([pd.read_csv(file, na_filter=False) for file in categories["b_factor_data"]]).to_csv("combined/b_factor_data.csv", index=False)
pd.concat([pd.read_csv(file, na_filter=False) for file in categories["don_hbonds_nr"]]).to_csv("combined/don_hbonds_nr.csv", index=False)
pd.concat([pd.read_csv(file, na_filter=False) for file in categories["prot_don_hbonds_nr"]]).to_csv("combined/prot_don_hbonds_nr.csv", index=False)
pd.concat([pd.read_csv(file, na_filter=False) for file in categories["acc_hbonds_nr"]]).to_csv("combined/acc_hbonds_nr.csv", index=False)
pd.concat([pd.read_csv(file, na_filter=False) for file in categories["hbond_data"]]).to_csv("combined/hbond_data.csv", index=False)
