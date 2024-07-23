"""
This script will eventually do something.
"""

import sys
import pandas as pd
import residue_library

# Redirect stdout and stderr to log files.
stdout = sys.stdout
stderr = sys.stderr
stdout_file = open(snakemake.log.stdout, mode='w')
stderr_file = open(snakemake.log.stderr, mode='w')
sys.stdout = stdout_file
sys.stderr = stderr_file

# Initialize the combined dataframe.
combined_df = pd.DataFrame(columns=['don_index', 'don_name', 'don_resn', 'don_resi', 'don_chain', 'don_segi', 'count_1',
                                    'count_2', 'b-factor', 'DOI', 'don_can_NA', 'acc_index', 'acc_name', 'acc_resn',
                                    'acc_resi', 'acc_chain', 'acc_segi', 'don_resn_name', 'acc_resn_name',
                                    'don_acc_distance', 'h_acc_distance', 'don_angle', 'h_angle', 'h_dihedral',
                                    'h_name', 'AOI', 'model', 'PDB', 'eq_class_member'])

# Combine the dataframes from the different equivalence class members. Do not incorporate empty dataframes.
for idx in range(len(snakemake.input.data)):
    try:
        member_df = pd.read_csv(snakemake.input.data[idx], comment="#", keep_default_na=False, na_values="NaN",
                                dtype={"don_resi": "object", "acc_resi": "object", "don_chain": "object",
                                       "acc_chain": "object"})
        if combined_df.empty and not member_df.empty:
            combined_df = member_df.copy()
        elif not member_df.empty:
            combined_df = pd.concat([combined_df, member_df])
    except (FileNotFoundError, pd.errors.EmptyDataError):
        continue
combined_df = combined_df.reset_index(drop=True)

# Prepare a list of acceptor residue and atom names from canonical amino acids and nucleic acids that have substantially
# greater negative charge.
can_neg = []
can_neut = []
for residue in residue_library.RESIDUE_LIBRARY:
    for acceptor in residue['acc']:
        if acceptor[3] < 0:
            can_neg.append(f"{residue['res']}.{acceptor[0]}")
        elif residue['res'] == 'ASP' and acceptor[0] == 'OD1':
            can_neg.append(f"{residue['res']}.{acceptor[0]}")
        elif residue['res'] == 'GLU' and acceptor[0] == 'OE1':
            can_neg.append(f"{residue['res']}.{acceptor[0]}")
        elif residue['res'] in ['A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT'] and acceptor[0] == 'OP1':
            can_neg.append(f"{residue['res']}.{acceptor[0]}")
        else:
            can_neut.append(f"{residue['res']}.{acceptor[0]}")

# Add a new column to identify acceptors from canonical amino acids and nucleic acids that have substantially greater
# negative charges.
combined_df["acc_charge"] = pd.NA
combined_df.loc[~combined_df.isna()['acc_index'], "acc_charge"] = "other"
combined_df.loc[combined_df["acc_resn_name"].isin(can_neg), "acc_charge"] = "can_neg"
combined_df.loc[combined_df["acc_resn_name"].isin(can_neut), "acc_charge"] = "can_neut"
combined_df = combined_df.drop(columns=["don_resn_name", "acc_resn_name"])

# Add a new column to identify donor-acceptor pairs that will be used for evaluating H-bonding geometry. Each atom pair
# involves a donor of interest.
combined_df = combined_df.assign(geom=0)
don_acc_grp = ["don_index", "acc_index", "eq_class_member"]
(combined_df.loc[
    # Only include atom pairs involving donors of interest.
    (combined_df["DOI"] == 1) &
    # For a given donor-acceptor pair, only include the hydrogen with the smaller D-H...A distance.
    (combined_df.groupby(don_acc_grp)["h_acc_distance"]
     .transform(lambda grp: [mem == grp.min() for mem in grp])) &
    # Do not consider A(N6)-U(O4), C(N4)-G(O6), G(N2)-C(O2), or G(N2)-C(N3) atom pairs.
    (~combined_df[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["A", "N6", "U", "O4"])
     .all(axis='columns')) &
    (~combined_df[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["C", "N4", "G", "O6"])
     .all(axis='columns')) &
    (~combined_df[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["G", "N2", "C", "O2"])
     .all(axis='columns')) &
    (~combined_df[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["G", "N2", "C", "N3"])
     .all(axis='columns')), "geom"]) = 1

# Set the H-bonding criteria.
H_DIST_MAX = snakemake.config["h_dist_max"]
H_ANG_MIN = snakemake.config["h_ang_min"]
DON_DIST_MAX = snakemake.config["don_dist_max"]
DON_ANG_MIN = snakemake.config["don_ang_min"]
DON_ANG_MAX = snakemake.config["don_ang_max"]


# Define a function to determine whether a donor/acceptor atom pair meets the H-bonding criteria.
def h_bonding(row):
    # If there are no H-bonding measurements, return NaN.
    if pd.isna(row["don_acc_distance"]) and pd.isna(row["h_acc_distance"]):
        return row
    # If the only distance measurement is between donor and acceptor, the donor is rotatable.
    elif pd.isna(row["h_acc_distance"]):
        if ((row["don_acc_distance"] <= DON_DIST_MAX) & (row["don_angle"] >= DON_ANG_MIN) &
                (row["don_angle"] <= DON_ANG_MAX)):
            row["h_bond"] = 1
            return row
        else:
            row["h_bond"] = 0
            return row
    # If both distance measurements exist, the donor is not rotatable.
    else:
        if (row["h_acc_distance"] <= H_DIST_MAX) & (row["h_angle"] >= H_ANG_MIN):
            row["h_bond"] = 1
            return row
        else:
            row["h_bond"] = 0
            return row


# Add a new column to identify whether each donor-acceptor pair meets the H-bond criteria.
combined_df = combined_df.assign(h_bond=pd.NA)
combined_df = combined_df.apply(h_bonding, axis=1)


# Define a function to classify donors of interest as no (0), single (1), or dual (2), H-bonding. The donors that meet
# the dual classification must use two of its hydrogens in the H-bonding interactions which must involve at least two
# different acceptors.
def classify(grp):
    h_name_list = []
    acc_index_list = []
    for row in grp.itertuples():
        if not pd.isna(row.h_bond):
            if row.h_bond == 1 and row.h_name not in h_name_list:
                h_name_list.append(row.h_name)
            if row.h_bond == 1 and row.acc_index not in acc_index_list:
                acc_index_list.append(row.acc_index)
    if len(h_name_list) == 1 or len(acc_index_list) == 1:
        grp['type'] = 1
        return grp
    elif len(h_name_list) == 2 and len(acc_index_list) >= 2:
        grp['type'] = 2
        return grp
    else:
        grp['type'] = 0
        return grp


# Add a new column that classifies the donors of interest as being no (0), single (1), or dual (2), H-bonding.
combined_df = combined_df.assign(type=pd.NA)
combined_df.loc[combined_df['DOI'] == 1, ~combined_df.columns.isin(["don_index", "eq_class_member"])] = (
    combined_df[combined_df['DOI'] == 1]
    .groupby(["don_index", "eq_class_member"], group_keys=False)
    .apply(classify, include_groups=False))

# Create dataframes of specific H-bonding interactions.
n6_o4_h_bond = (combined_df[
                    (combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                     .eq(["N6", "A", "O4", "U", 1]).all(axis='columns')) |
                    (combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                     .eq(["N6", "DA", "O4", "U", 1]).all(axis='columns')) |
                    (combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                     .eq(["N6", "A", "O4", "DT", 1]).all(axis='columns')) |
                    (combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                     .eq(["N6", "DA", "O4", "DT", 1]).all(axis='columns'))
                ].rename(columns=lambda col: f'{col}_N6_O4'))
n3_n1_h_bond = (combined_df[
                    (combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                     .eq(["N3", "U", "N1", "A", 1]).all(axis='columns')) |
                    (combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                     .eq(["N3", "DT", "N1", "A", 1]).all(axis='columns')) |
                    (combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                     .eq(["N3", "U", "N1", "DA", 1]).all(axis='columns')) |
                    (combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                     .eq(["N3", "DT", "N1", "DA", 1]).all(axis='columns'))
                ].rename(columns=lambda col: f'{col}_N3_N1'))
n4_o6_h_bond = (combined_df[
                    (combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                     .eq(["N4", "C", "O6", "G", 1]).all(axis='columns')) |
                    (combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                     .eq(["N4", "DC", "O6", "G", 1]).all(axis='columns')) |
                    (combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                     .eq(["N4", "C", "O6", "DG", 1]).all(axis='columns')) |
                    (combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                     .eq(["N4", "DC", "O6", "DG", 1]).all(axis='columns'))
                ].rename(columns=lambda col: f'{col}_N4_O6'))
n1_n3_h_bond = (combined_df[
                    (combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                     .eq(["N1", "G", "N3", "C", 1]).all(axis='columns')) |
                    (combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                     .eq(["N1", "DG", "N3", "C", 1]).all(axis='columns')) |
                    (combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                     .eq(["N1", "G", "N3", "DC", 1]).all(axis='columns')) |
                    (combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                     .eq(["N1", "DG", "N3", "DC", 1]).all(axis='columns'))
                ].rename(columns=lambda col: f'{col}_N1_N3'))
n2_o2_h_bond = (combined_df[
                    (combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                     .eq(["N2", "G", "O2", "C", 1]).all(axis='columns')) |
                    (combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                     .eq(["N2", "DG", "O2", "C", 1]).all(axis='columns')) |
                    (combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                     .eq(["N2", "G", "O2", "DC", 1]).all(axis='columns')) |
                    (combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                     .eq(["N2", "DG", "O2", "DC", 1]).all(axis='columns'))
                ].rename(columns=lambda col: f'{col}_N2_O2'))

# Prepare a dataframe of nucleobases containing a donor of interest and involved in a canonical base pair.
# AU base pair
base_pair_df = (n6_o4_h_bond[n6_o4_h_bond["DOI_N6_O4"] == 1]
                .merge(n3_n1_h_bond,
                       left_on=["don_resn_N6_O4", "don_resi_N6_O4", "don_chain_N6_O4", "acc_resn_N6_O4",
                                "acc_resi_N6_O4", "acc_chain_N6_O4"],
                       right_on=["acc_resn_N3_N1", "acc_resi_N3_N1", "acc_chain_N3_N1", "don_resn_N3_N1",
                                 "don_resi_N3_N1", "don_chain_N3_N1"], how='inner')
                .rename(columns={"don_index_N6_O4": "don_index", "eq_class_member_N6_O4": "eq_class_member"})
                .assign(base_pair="AU").loc[:, ["don_index", "eq_class_member", "base_pair"]])
# CG base pair
base_pair_df = pd.concat([base_pair_df,
                          (n4_o6_h_bond[n4_o6_h_bond["DOI_N4_O6"] == 1]
                           .merge(n1_n3_h_bond,
                                  left_on=["don_resn_N4_O6", "don_resi_N4_O6", "don_chain_N4_O6", "acc_resn_N4_O6",
                                           "acc_resi_N4_O6", "acc_chain_N4_O6"],
                                  right_on=["acc_resn_N1_N3", "acc_resi_N1_N3", "acc_chain_N1_N3", "don_resn_N1_N3",
                                            "don_resi_N1_N3", "don_chain_N1_N3"], how='inner')
                           .merge(n2_o2_h_bond,
                                  left_on=["don_resn_N4_O6", "don_resi_N4_O6", "don_chain_N4_O6", "acc_resn_N4_O6",
                                           "acc_resi_N4_O6", "acc_chain_N4_O6"],
                                  right_on=["acc_resn_N2_O2", "acc_resi_N2_O2", "acc_chain_N2_O2", "don_resn_N2_O2",
                                            "don_resi_N2_O2", "don_chain_N2_O2"], how='inner')
                           .rename(columns={"don_index_N4_O6": "don_index", "eq_class_member_N4_O6": "eq_class_member"})
                           .assign(base_pair="CG").loc[:, ["don_index", "eq_class_member", "base_pair"]])])
# GC base pair
base_pair_df = pd.concat([base_pair_df,
                          (n2_o2_h_bond[n2_o2_h_bond["DOI_N2_O2"] == 1]
                           .merge(n1_n3_h_bond,
                                  left_on=["don_resn_N2_O2", "don_resi_N2_O2", "don_chain_N2_O2", "acc_resn_N2_O2",
                                           "acc_resi_N2_O2", "acc_chain_N2_O2"],
                                  right_on=["don_resn_N1_N3", "don_resi_N1_N3", "don_chain_N1_N3", "acc_resn_N1_N3",
                                            "acc_resi_N1_N3", "acc_chain_N1_N3"], how='inner')
                           .merge(n4_o6_h_bond,
                                  left_on=["don_resn_N2_O2", "don_resi_N2_O2", "don_chain_N2_O2", "acc_resn_N2_O2",
                                           "acc_resi_N2_O2", "acc_chain_N2_O2"],
                                  right_on=["acc_resn_N4_O6", "acc_resi_N4_O6", "acc_chain_N4_O6", "don_resn_N4_O6",
                                            "don_resi_N4_O6", "don_chain_N4_O6"], how='inner')
                           .rename(columns={"don_index_N2_O2": "don_index", "eq_class_member_N2_O2": "eq_class_member"})
                           .assign(base_pair="GC").loc[:, ["don_index", "eq_class_member", "base_pair"]])])
combined_df = combined_df.merge(base_pair_df, how='outer')

# Write the processed and combined data to a csv file.
combined_df.to_csv(snakemake.output.combined, index=False, na_rep='NaN')


# Define a function to determine whether the residue donating to the acceptor of interest is the same as the residue
# accepting from the donor of interest. If it is the same residue, the interactions are considered to be overlapping.
# The acceptor of interest and donor of interest are a part of the same nucleobase. Additionally, identify whether any
# of the atoms accepting an H-bond from the donor of interest bear substantial negative charge.
def check_overlap_charge(grp):
    overlap = 0
    neg_count = 0
    for row in grp.itertuples():
        # Is the residue accepting the H-bond from the amine the same residue that is donating an H-bond to the acceptor
        # of interest?
        if row.acc_resn == row.don_resn_AOI and row.acc_resi == row.don_resi_AOI and row.acc_chain == row.don_chain_AOI:
            # Are the acceptor and donor atoms both part of the nucleobase (e.i., not backbone atoms)?
            if row.acc_name not in ["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"] and row.don_name_AOI not in ["O2'"]:
                overlap = 1
        if row.acc_charge == "can_neg":
            neg_count += 1
    grp['overlap'] = overlap
    # The DOI does not donate to any anionic acceptors.
    if neg_count == 0:
        grp["neg_charge"] = "None"
    # The DOI donates to neutral and anionic acceptors.
    elif neg_count < len(grp):
        grp["neg_charge"] = "Partial"
    # The DOI only donates to anionic acceptors.
    elif neg_count == len(grp):
        grp["neg_charge"] = "All"
    return grp


# Define a function to determine whether and donor of interest donates to any residues not included in the
# included_residues list specified within the config file. If it does, exclude the donor of interest from further
# consideration by returning an empty data frame.
def check_included_res(grp):
    for row in grp.itertuples():
        if row.acc_resn not in snakemake.config["included_residues"]:
            grp = grp.truncate(after=-1)
            return grp
    return grp


# Create a data frame that includes H-bonding atom pairs that include an acceptor of interest.
aoi_df = (combined_df[(combined_df["AOI"] == 1) & (combined_df["h_bond"] == 1)]
          .drop_duplicates(subset=["acc_index", "don_index", "eq_class_member"])
          .loc[:, ['don_index', 'don_name', 'don_resn', 'don_resi', 'don_chain', 'acc_index', 'acc_name', 'acc_resn',
                   'acc_resi', 'acc_chain', 'eq_class_member', 'don_acc_distance']])

# Create a data frame that includes H-bonding atom pairs that include a donor of interest which only donates to residues
# included in the included_residues list specified in the config file.
doi_df = combined_df[(combined_df["DOI"] == 1) & (combined_df["h_bond"] == 1)]
doi_df.loc[:, ~doi_df.columns.isin(["don_index", "eq_class_member"])] = (
    doi_df.groupby(['don_index', 'eq_class_member'], group_keys=False).apply(check_included_res, include_groups=False))

# Create a merged data frame where residues bearing donors of interest are matched with residues bearing acceptors of
# interest.
merged_df = (doi_df.merge(aoi_df, left_on=['don_resn', 'don_resi', 'don_chain', 'eq_class_member'],
                          right_on=['acc_resn', 'acc_resi', 'acc_chain', 'eq_class_member'], suffixes=[None, "_AOI"]))

# Add overlap and neg_charge columns to the merged data frame. Refer to the check_overlap_charge function definition
# above for what these columns indicate.
merged_df.loc[:, ["overlap", "neg_charge"]] = pd.NA
merged_df.loc[:, ~merged_df.columns.isin(['don_index_AOI', "acc_index_AOI", "eq_class_member"])] = (
    merged_df.groupby(['don_index_AOI', 'acc_index_AOI', 'eq_class_member'], group_keys=False)
    .apply(check_overlap_charge, include_groups=False))

# Remove unnecessary columns and redundant rows from the merged data frame.
drop_columns = ["don_index", "don_resi", "don_chain", "don_segi", "count_1", "count_2", "b-factor", "DOI", "don_can_NA",
                "acc_index", "acc_name", "acc_resn", "acc_resi", "acc_chain", "acc_segi", "don_acc_distance",
                "h_acc_distance", "don_angle", "h_angle", "h_dihedral", "h_name", "AOI", "h_dihedral", "model", "PDB",
                "eq_class_member", "acc_charge", "geom", "h_bond", "don_index_AOI", "don_name_AOI", "don_resn_AOI",
                "don_resi_AOI", "don_chain_AOI", "acc_index_AOI", "acc_resn_AOI", "acc_resi_AOI", "acc_chain_AOI"]
merged_df = (merged_df.drop_duplicates(subset=['don_index_AOI', 'acc_index_AOI', 'eq_class_member'])
             .drop(columns=drop_columns))

# Write the merged data frame to a csv file.
merged_df.to_csv(snakemake.output.distances, index=False, na_rep='NaN')

# Close files and reset stdout and stderr.
stdout_file.close()
stderr_file.close()
sys.stdout = stdout
sys.stderr = stderr
