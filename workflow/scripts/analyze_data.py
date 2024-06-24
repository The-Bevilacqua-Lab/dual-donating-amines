"""
This script will eventually do something.
"""

import numpy as np
import pandas as pd

# Set the H-bonding criteria.
H_DIST_MAX = snakemake.config["h_dist_max"]
H_ANG_MIN = snakemake.config["h_ang_min"]

# Prepare a list of acceptor residue names and atom names that have substantially greater negative charge.
neg_acc_resn = []
neg_acc_name = []
for residue in const.RESIDUE_LIBRARY:
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

# Get the data.
count_data = pd.read_csv(snakemake.input.counts, keep_default_na=False, na_values="NaN", dtype={"resi": "object"})
h_bond_data = pd.read_csv(snakemake.input.h_bond, keep_default_na=False, na_values="NaN",
                          dtype={"don_resi": "object", "acc_resi": "object"})
nuc_data_raw = pd.read_csv(snakemake.input.nuc, keep_default_na=False, na_values="NaN", dtype={"resi": "object"})
b_factor_data = pd.read_csv(snakemake.input.b_factor, keep_default_na=False, na_values="NaN", dtype={"resi": "object"})

# Create a new dataframe of nucleobases containing a donor of interest where the mean of the b-factors are below the
# cutoff.
nuc_data = (nuc_data_raw.merge(b_factor_data[(b_factor_data["subset"] == "sidechain") &
                                             (b_factor_data["mean"] < snakemake.config["b_factor_cutoff"])],
                               on=["resn", "resi", "chain", "eq_class_members"], how='inner')
            .drop(columns=["mean", "subset"]))

# Prepare a list containing the residue and atom names of the donors of interest.
donors_of_interest = []
for donor in snakemake.config["donors_of_interest"]:
    donors_of_interest.append([donor.split(".")[0], donor.split(".")[1]])

# Identify atom pairs that meet the H-bond criteria and include a donor of interest. The nucleobase containing the donor
# of interest must also meet the b-factor criteria.
don_h_bonds = (h_bond_data[(h_bond_data["h_acc_distance"] <= H_DIST_MAX) &
                           (h_bond_data["h_angle"] >= H_ANG_MIN)]
               .merge(pd.DataFrame(donors_of_interest, columns=["don_resn", "don_name"]), how='inner')
               .merge(nuc_data, left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                      right_on=["resn", "resi", "chain", "eq_class_members"], how='inner')
               .drop(columns=["resn", "resi", "chain"]))

# Write data on donor-acceptor pairs involving donors of interest to csv files.
don_acc_grp = ["don_index", "acc_index", "eq_class_members"]
don_atom_pairs = (don_h_bonds[
     # Only include hydrogens with smaller D-H...A distances.
     (don_h_bonds.groupby(don_acc_grp)["h_acc_distance"]
      .transform(lambda grp: [mem == grp.min() for mem in grp])) &
     # Do not consider A(N6)-U(O4), C(N4)-G(O6), G(N2)-C(O2), or G(N2)-C(N3) atom pairs.
     (~don_h_bonds[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["A", "N6", "U", "O4"])
      .all(axis='columns')) &
     (~don_h_bonds[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["C", "N4", "G", "O6"])
      .all(axis='columns')) &
     (~don_h_bonds[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["G", "N2", "C", "O2"])
      .all(axis='columns')) &
     (~don_h_bonds[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["G", "N2", "C", "N3"])
      .all(axis='columns'))])
don_atom_pairs.to_csv(snakemake.output.atom_pairs, index=False, columns=["don_acc_distance", "h_acc_distance", "h_angle"])

# Identify nucleobases that donate either one or at least two H-bonds via their exocyclic amines. The H-bonds from the
# latter nucleobases must involve both exocyclic amine hydrogens and at least two different acceptors.
single_don_h_bonds = don_h_bonds[
    (don_h_bonds.groupby(["don_index", "eq_class_members"])["h_name"].transform("nunique") == 1) |
    (don_h_bonds.groupby(["don_index", "eq_class_members"])["acc_index"].transform("nunique") == 1)]
dual_don_h_bonds = don_h_bonds[
    (don_h_bonds.groupby(["don_index", "eq_class_members"])["h_name"].transform("nunique") == 2) &
    (don_h_bonds.groupby(["don_index", "eq_class_members"])["acc_index"].transform("nunique") >= 2)]

# Prepare dataframes of no, single, and dual H-bonding donors of interest with just residue information.
single_don_res = (single_don_h_bonds.loc[:, ["eq_class_members", "don_resn", "don_resi", "don_chain"]]
                  .drop_duplicates().rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"}))
dual_don_res = (dual_don_h_bonds.loc[:, ["eq_class_members", "don_resn", "don_resi", "don_chain"]]
                .drop_duplicates().rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"}))
no_don_res = pd.concat([single_don_res, dual_don_res, nuc_data]).drop_duplicates(keep=False)

# Prepare a dataframe of residues that contain donors of interest with nearby heavy atom counts and H-bonding type
# (no, single, or dual) information.
no_count = no_don_res.merge(count_data, how='inner').drop(columns=["index", "name"])
single_count = single_don_res.merge(count_data, how='inner').drop(columns=["index", "name"])
dual_count = dual_don_res.merge(count_data, how='inner').drop(columns=["index", "name"])
no_count["type"] = 0
single_count["type"] = 1
dual_count["type"] = 2
organized_count = pd.concat([no_count, single_count, dual_count])
organized_count["volume_1"] = (4/3)*np.pi*np.power(snakemake.config["count_dist_1"], 3)
organized_count["volume_2"] = (4/3)*np.pi*np.power(snakemake.config["count_dist_2"], 3)
organized_count.to_csv(snakemake.output.counts, index=False)
