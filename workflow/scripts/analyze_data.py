"""
This script will eventually do something.
"""

import numpy as np
import pandas as pd
import residue_library

# Initialize the combined dataframe.
combined_df = pd.DataFrame(columns=["don_index", "don_name", "don_resn", "don_resi", "don_chain", "don_segi", "count_1",
                                    "count_2", "b-factor", "DOI", "acc_index", "acc_name", "acc_resn", "acc_resi",
                                    "acc_chain", "acc_segi", "don_acc_distance", "h_acc_distance", "h_angle",
                                    "h_dihedral", "h_name", "model", "PDB", "eq_class_member"])

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
    except FileNotFoundError:
        continue
combined_df = combined_df.reset_index(drop=True)

# Filter out nucleobases containing donors of interest that do not meet the b-factor criteria.
combined_df = combined_df[(combined_df["DOI"] == 0) | ((combined_df["DOI"] == 1) &
                                                       (combined_df["b-factor"] < snakemake.config["b_factor_cutoff"]))]

# Write data on donor-acceptor pairs involving donors of interest to csv files.
don_acc_grp = ["don_index", "acc_index", "eq_class_member"]
(combined_df[
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
     .all(axis='columns'))]
 .to_csv(snakemake.output.atom_pairs, index=False, columns=["don_acc_distance", "h_acc_distance", "h_angle"]))

# Set the H-bonding criteria.
H_DIST_MAX = snakemake.config["h_dist_max"]
H_ANG_MIN = snakemake.config["h_ang_min"]

# Add a new column to identify whether each donor-acceptor pair meets the H-bond criteria.
combined_df["h_bond"] = (combined_df["h_acc_distance"] <= H_DIST_MAX) & (combined_df["h_angle"] >= H_ANG_MIN)


# Define a function to classify donors of interest as no (0), single (1), or dual (2), H-bonding. The donors that meet
# the dual classification must use two of its hydrogens in the H-bonding interactions which must involve at least two
# different acceptors.
def classify(grp):
    h_name_list = []
    acc_index_list = []
    for row in grp.itertuples():
        if not pd.isna(row.h_bond):
            if row.h_bond and row.h_name not in h_name_list:
                h_name_list.append(row.h_name)
            if row.h_bond and row.acc_index not in acc_index_list:
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
combined_df.loc[combined_df['DOI'] == 1, combined_df.columns != "don_index"] = (combined_df[combined_df['DOI'] == 1]
                                                                                .groupby('don_index', group_keys=False)
                                                                                .apply(classify, include_groups=False))

# Write data on residues that contain donors of interest with nearby heavy atom counts and H-bonding type (no, single,
# or dual) information.
(combined_df.loc[combined_df['DOI'] == 1, ["don_index", "don_resn", "count_1", "count_2", "eq_class_member", "type"]]
 .drop_duplicates(subset=["don_index", "eq_class_member"])
 .assign(volume_1=(4/3)*np.pi*snakemake.config["count_dist_1"]**3,
         volume_2=(4/3)*np.pi*snakemake.config["count_dist_2"]**3).to_csv(snakemake.output.counts, index=False))

# Create dataframes of specific H-bonding interactions.
n6_o4_h_bond = (combined_df[combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                .eq(["N6", "A", "O4", "U", True]).all(axis='columns')].rename(columns=lambda col: f'{col}_N6_O4'))
n3_n1_h_bond = (combined_df[combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                .eq(["N3", "U", "N1", "A", True]).all(axis='columns')].rename(columns=lambda col: f'{col}_N3_N1'))
n4_o6_h_bond = (combined_df[combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                .eq(["N4", "C", "O6", "G", True]).all(axis='columns')].rename(columns=lambda col: f'{col}_N4_O6'))
n1_n3_h_bond = (combined_df[combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                .eq(["N1", "G", "N3", "C", True]).all(axis='columns')].rename(columns=lambda col: f'{col}_N1_N3'))
n2_o2_h_bond = (combined_df[combined_df[["don_name", "don_resn", "acc_name", "acc_resn", "h_bond"]]
                .eq(["N2", "G", "O2", "C", True]).all(axis='columns')].rename(columns=lambda col: f'{col}_N2_O2'))

# Write a csv of nucleobases containing a donor of interest and involved in a canonical base pair.
# AU base pair
(n6_o4_h_bond[n6_o4_h_bond["DOI_N6_O4"] == 1]
 .merge(n3_n1_h_bond,
        left_on=["don_resn_N6_O4", "don_resi_N6_O4", "don_chain_N6_O4", "acc_resn_N6_O4", "acc_resi_N6_O4",
                 "acc_chain_N6_O4"],
        right_on=["acc_resn_N3_N1", "acc_resi_N3_N1", "acc_chain_N3_N1", "don_resn_N3_N1", "don_resi_N3_N1",
                  "don_chain_N3_N1"], how='inner').assign(base_pair="AU")
 .to_csv(snakemake.output.base_pair, index=False))
# CG base pair
(n4_o6_h_bond[n4_o6_h_bond["DOI_N4_O6"] == 1]
 .merge(n1_n3_h_bond,
        left_on=["don_resn_N4_O6", "don_resi_N4_O6", "don_chain_N4_O6", "acc_resn_N4_O6", "acc_resi_N4_O6",
                 "acc_chain_N4_O6"],
        right_on=["acc_resn_N1_N3", "acc_resi_N1_N3", "acc_chain_N1_N3", "don_resn_N1_N3", "don_resi_N1_N3",
                  "don_chain_N1_N3"], how='inner')
 .merge(n2_o2_h_bond,
        left_on=["don_resn_N4_O6", "don_resi_N4_O6", "don_chain_N4_O6", "acc_resn_N4_O6", "acc_resi_N4_O6",
                 "acc_chain_N4_O6"],
        right_on=["acc_resn_N2_O2", "acc_resi_N2_O2", "acc_chain_N2_O2", "don_resn_N2_O2", "don_resi_N2_O2",
                  "don_chain_N2_O2"], how='inner').assign(base_pair="CG")
 .to_csv(snakemake.output.base_pair, index=False, mode='a'))
# GC base pair
(n2_o2_h_bond[n2_o2_h_bond["DOI_N2_O2"] == 1]
 .merge(n1_n3_h_bond,
        left_on=["don_resn_N2_O2", "don_resi_N2_O2", "don_chain_N2_O2", "acc_resn_N2_O2", "acc_resi_N2_O2",
                 "acc_chain_N2_O2"],
        right_on=["don_resn_N1_N3", "don_resi_N1_N3", "don_chain_N1_N3", "acc_resn_N1_N3", "acc_resi_N1_N3",
                  "acc_chain_N1_N3"], how='inner')
 .merge(n4_o6_h_bond,
        left_on=["don_resn_N2_O2", "don_resi_N2_O2", "don_chain_N2_O2", "acc_resn_N2_O2", "acc_resi_N2_O2",
                 "acc_chain_N2_O2"],
        right_on=["acc_resn_N4_O6", "acc_resi_N4_O6", "acc_chain_N4_O6", "don_resn_N4_O6", "don_resi_N4_O6",
                  "don_chain_N4_O6"], how='inner').assign(base_pair="GC")
 .to_csv(snakemake.output.base_pair, index=False, mode='a'))

# Prepare a list of acceptor residue names and atom names that have substantially greater negative charge.
neg_acc_resn = []
neg_acc_name = []
for residue in residue_library.RESIDUE_LIBRARY:
    for acceptor in residue['acc']:
        if acceptor[3] < 0:
            neg_acc_resn.append(residue['res'])
            neg_acc_name.append(acceptor[0])
        elif residue['res'] == 'ASP' and acceptor[0] == 'OD1':
            neg_acc_resn.append(residue['res'])
            neg_acc_name.append(acceptor[0])
        elif residue['res'] == 'GLU' and acceptor[0] == 'OE1':
            neg_acc_resn.append(residue['res'])
            neg_acc_name.append(acceptor[0])
        elif residue['res'] in ['A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT'] and acceptor[0] == 'OP1':
            neg_acc_resn.append(residue['res'])
            neg_acc_name.append(acceptor[0])

# Prepare a dataframe of acceptors that have substantially greater negative charges.
neg_acceptors = pd.DataFrame({
    "acc_resn": neg_acc_resn,
    "acc_name": neg_acc_name
})
