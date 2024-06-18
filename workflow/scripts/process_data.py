"""
This module will eventually do something.
"""

import pandas as pd
import residue_library

# Set the H-bonding criteria.
H_DIST_MAX = snakemake.config["h_dist_max"]
H_ANG_MIN = snakemake.config["h_ang_min"]

# Construct a list of tuples containing atom and residue names that describe atoms capable of both donating and
# accepting an H-bond.
don_acc_atoms = []
for residue in const.RESIDUE_LIBRARY:
    for donor in residue['don']:
        for acceptor in residue['acc']:
            if donor[0] == acceptor[0]:
                don_acc_atoms.append((residue['res'], donor[0]))


# Organize the data and apply the H-bonding criteria.
def process_data(h_bond_file, nuc_file, b_factor_file, eq_class_members_file):

    # Extract the data from the h_bond csv file and remove redundant lines.
    unique_col_comb = ["don_index", "acc_index", "h_name"]
    # noinspection PyTypeChecker
    h_bond_data = (pd.read_csv(h_bond_file, comment="#", keep_default_na=False,
                               na_values={"h_acc_distance": "NaN", "h_angle": "NaN", "h_dihedral": "NaN"},
                               dtype={"don_chain": "object", "acc_chain": "object"})
                   .drop_duplicates(subset=unique_col_comb))

    # Extract the data from the nuc csv file.
    nuc_data = pd.read_csv(nuc_file, comment="#", na_filter=False, dtype={"chain": "object"})

    # Create a dictionary based on the file containing information on nucleobase b-factors.
    b_factor_grp = ["resn", "resi", "chain", "subset"]
    b_factor_data = pd.read_csv(b_factor_file, comment="#", na_filter=False, dtype={"chain": "object"})
    b_factor_data["mean"] = b_factor_data.groupby(b_factor_grp)["b-factor"].transform("mean")
    b_factor_data = b_factor_data.drop_duplicates(subset=b_factor_grp).drop(columns=["index", "name", "b-factor"])

    # Extract the data from the equivalence class file.
    eq_class_members_data = (pd.read_csv(eq_class_members_file, header=None, comment="#",
                                         # TODO "True" should eventually be replaced with "snakemake.config["commit_hash"]"
                                         skiprows=[3 if True else 2],
                                         na_filter=False, dtype="object")
                             .set_index(pd.Series(["pdb_id", "model", "chain"])).transpose())

    # Identify atom pairs that meet the H-bond criteria and include a donor of interest.
    don_h_bonds = (h_bond_data[(h_bond_data["h_acc_distance"] <= H_DIST_MAX) &
                               (h_bond_data["h_angle"] >= H_ANG_MIN)]
                   .merge(pd.DataFrame(const.DONORS_OF_INTEREST, columns=["don_resn", "don_name"]), how='inner')
                   .merge(eq_class_members_data, left_on="don_chain", right_on="chain", how='inner')
                   .drop(columns=["pdb_id", "model", "chain"]))

    # Identify atom pairs that meet the H-bond criteria and include a protonated donor of interest. Values in the _merge
    # column matching left_only indicate acceptors that cannot typically also donate an H-bond.
    prot_don_h_bonds_unfiltered = (h_bond_data[(h_bond_data["h_acc_distance"] <= H_DIST_MAX) &
                                               (h_bond_data["h_angle"] >= H_ANG_MIN)]
                                   .merge(pd.DataFrame(const.PROT_DONORS_OF_INTEREST, columns=["don_resn", "don_name"]),
                                          how='inner')
                                   .merge(pd.DataFrame(don_acc_atoms, columns=["acc_resn", "acc_name"]), how='left',
                                          indicator=True))

    # Filter prot_don_h_bonds_unfiltered such that only acceptors that cannot typically also donate an H-bond are
    # included.
    prot_don_h_bonds = (prot_don_h_bonds_unfiltered[prot_don_h_bonds_unfiltered["_merge"] == "left_only"]
                        .drop(columns="_merge")
                        .merge(eq_class_members_data, left_on="don_chain", right_on="chain", how='inner')
                        .drop(columns=["pdb_id", "model", "chain"]))

    # Identify atom pairs that meet the H-bond criteria and include an acceptor of interest.
    acc_h_bonds = (h_bond_data[(h_bond_data["h_acc_distance"].notna()) & (h_bond_data["h_acc_distance"] <= H_DIST_MAX) &
                               (h_bond_data["h_angle"] >= H_ANG_MIN)]
                   .merge(pd.DataFrame(const.ACCEPTORS_OF_INTEREST, columns=["acc_resn", "acc_name"]), how='inner')
                   .merge(eq_class_members_data, left_on="acc_chain", right_on="chain", how='inner')
                   .drop(columns=["pdb_id", "model", "chain"]))

    # Specify the ID of the equivalence class member within an additional column in the dataframes.
    eq_class_member_id = h_bond_file[13:-11]
    h_bond_data["eq_class_members"] = eq_class_member_id
    nuc_data["eq_class_members"] = eq_class_member_id
    b_factor_data["eq_class_members"] = eq_class_member_id
    don_h_bonds["eq_class_members"] = eq_class_member_id
    prot_don_h_bonds["eq_class_members"] = eq_class_member_id
    acc_h_bonds["eq_class_members"] = eq_class_member_id

    # Return the dataframes.
    return h_bond_data, nuc_data, b_factor_data, don_h_bonds, prot_don_h_bonds, acc_h_bonds


# Initialize combined dataframes.
h_bond_combined = pd.DataFrame()
nuc_combined = pd.DataFrame()
b_factor_combined = pd.DataFrame()
don_h_bonds_combined = pd.DataFrame()
prot_don_h_bonds_combined = pd.DataFrame()
acc_h_bonds_combined = pd.DataFrame()

# Run process_data() and combine the resulting dataframes from different equivalence class members.
for idx in range(len(snakemake.input.h_bond)):
    if idx == 0:
        h_bond_df, nuc_df, b_factor_df, don_h_bonds_df, prot_don_h_bonds_df, acc_h_bonds_df = process_data(
            snakemake.input.h_bond[idx], snakemake.input.nuc[idx], snakemake.input.b_factor[idx],
            snakemake.input.eq_class_members[idx])
        h_bond_combined = h_bond_df.copy()
        nuc_combined = nuc_df.copy()
        b_factor_combined = b_factor_df.copy()
        don_h_bonds_combined = don_h_bonds_df.copy()
        prot_don_h_bonds_combined = prot_don_h_bonds_df.copy()
        acc_h_bonds_combined = acc_h_bonds_df.copy()
    else:
        h_bond_df, nuc_df, b_factor_df, don_h_bonds_df, prot_don_h_bonds_df, acc_h_bonds_df = process_data(
            snakemake.input.h_bond[idx], snakemake.input.nuc[idx], snakemake.input.b_factor[idx],
            snakemake.input.eq_class_members[idx])
        h_bond_combined = pd.concat([h_bond_combined, h_bond_df], ignore_index=True)
        nuc_combined = pd.concat([nuc_combined, nuc_df], ignore_index=True)
        b_factor_combined = pd.concat([b_factor_combined, b_factor_df], ignore_index=True)
        don_h_bonds_combined = pd.concat([don_h_bonds_combined, don_h_bonds_df], ignore_index=True)
        prot_don_h_bonds_combined = pd.concat([prot_don_h_bonds_combined, prot_don_h_bonds_df], ignore_index=True)
        acc_h_bonds_combined = pd.concat([acc_h_bonds_combined, acc_h_bonds_df], ignore_index=True)

# Write to csv files.
h_bond_combined.to_csv(snakemake.output.h_bond, index=False, na_rep='NaN')
nuc_combined.to_csv(snakemake.output.nuc, index=False, na_rep='NaN')
b_factor_combined.to_csv(snakemake.output.b_factor, index=False, na_rep='NaN')
don_h_bonds_combined.to_csv(snakemake.output.don_h_bonds, index=False, na_rep='NaN')
prot_don_h_bonds_combined.to_csv(snakemake.output.prot_don_h_bonds, index=False, na_rep='NaN')
acc_h_bonds_combined.to_csv(snakemake.output.acc_h_bonds, index=False, na_rep='NaN')
