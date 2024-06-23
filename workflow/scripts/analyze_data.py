"""
This script will eventually do something.
"""

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
# TODO prot donor should not donate to the same atom the exo amine is donating to (e.g., C788 A5 in 6xu8)
# noinspection PyTypeChecker
prot_don_h_bonds = pd.read_csv(snakemake.input.prot_don_h_bonds, na_filter=False,
                               dtype={"don_resi": "object", "acc_resi": "object"})
# noinspection PyTypeChecker
acc_h_bonds = pd.read_csv(snakemake.input.acc_h_bonds, keep_default_na=False,
                          na_values={"h_acc_distance": "NaN", "h_angle": "NaN", "h_dihedral": "NaN"},
                          dtype={"don_resi": "object", "acc_resi": "object"})
nuc_data_raw = pd.read_csv(snakemake.input.nuc, na_filter=False, dtype={"resi": "object"})
b_factor_data = pd.read_csv(snakemake.input.b_factor, na_filter=False, dtype={"resi": "object"})
# noinspection PyTypeChecker
h_bond_data = pd.read_csv(snakemake.input.h_bond, keep_default_na=False,
                          na_values={"h_acc_distance": "NaN", "h_angle": "NaN", "h_dihedral": "NaN"},
                          dtype={"don_resi": "object", "acc_resi": "object"})

# Create a new dataframe of nucleobases where the mean of the b-factors are below 79.
nuc_data = ((nuc_data_raw.merge(b_factor_data[(b_factor_data["subset"] == "sidechain") &
                                              (b_factor_data["mean"] < snakemake.config["b_factor_cutoff"])],
                                on=["resn", "resi", "chain", "eq_class_members"], how='inner').drop(columns=["mean"]))
            .drop(columns=["subset"]))

# Write data on donor-acceptor pairs to csv files.
don_acc_grp = ["don_index", "acc_index", "eq_class_members"]
a_n6_don = (h_bond_data[
     # Only include hydrogens with smaller D-H...A distances.
     (h_bond_data.groupby(don_acc_grp)["h_acc_distance"]
      .transform(lambda grp: [mem == grp.min() for mem in grp])) &
     # Do not consider A(N6)-U(O4) atom pairs.
     (~h_bond_data[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["A", "N6", "U", "O4"])
      .all(axis='columns')) &
     # Only consider atom pairs where A(N6) is the donor.
     (h_bond_data[["don_resn", "don_name"]].eq(["A", "N6"])
      .all(axis='columns'))])
a_n6_don["atom"] = "A(N6) Donor"
a_n1_acc = (h_bond_data[
     # Only include hydrogens with smaller D-H...A distances.
     (h_bond_data.groupby(don_acc_grp)["h_acc_distance"]
      .transform(lambda grp: [mem == grp.min() for mem in grp])) &
     # Only consider non-rotatable donors.
     (h_bond_data["h_acc_distance"].notna()) &
     # Do not consider U(N3)-A(N1) atom pairs.
     (~h_bond_data[["don_resn", "don_name", "acc_resn", "acc_name"]].eq(["U", "N3", "A", "N1"])
      .all(axis='columns')) &
     # Only consider atom pairs where A(N1) is the acceptor.
     (h_bond_data[["acc_resn", "acc_name"]].eq(["A", "N1"])
      .all(axis='columns'))])
a_n1_acc["atom"] = "A(N1) Acceptor"
pd.concat([a_n6_don, a_n1_acc], ignore_index=True).to_csv(snakemake.output.heatmap, index=False,
                                                          columns=["don_acc_distance", "don_angle", "acc_angle",
                                                                   "h_acc_distance", "h_angle", "atom"])

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

# Identify nucleobases that donate either one or at least two H-bonds via their exocyclic amines. The H-bonds from the
# latter nucleobases must involve both exocyclic amine hydrogens and at least two different acceptors.
single_don_h_bonds = don_h_bonds[
    (don_h_bonds.groupby(["don_index", "eq_class_members"])["h_name"].transform("nunique") == 1) |
    (don_h_bonds.groupby(["don_index", "eq_class_members"])["acc_index"].transform("nunique") == 1)]
dual_don_h_bonds = don_h_bonds[
    (don_h_bonds.groupby(["don_index", "eq_class_members"])["h_name"].transform("nunique") == 2) &
    (don_h_bonds.groupby(["don_index", "eq_class_members"])["acc_index"].transform("nunique") >= 2)]

# Prepare dataframes of single and dual H-bonding A(N6) with H-bond information.
a_n6_single_h_bond = (single_don_h_bonds[single_don_h_bonds["don_resn"] == "A"]
                      .merge(nuc_data,
                             left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                             right_on=["resn", "resi", "chain", "eq_class_members"],
                             how='inner').drop(columns=["resn", "resi", "chain"]))
a_n6_dual_h_bond = (dual_don_h_bonds[dual_don_h_bonds["don_resn"] == "A"]
                    .merge(nuc_data,
                           left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                           right_on=["resn", "resi", "chain", "eq_class_members"],
                           how='inner').drop(columns=["resn", "resi", "chain"]))

# Prepare dataframes of no, single, and dual H-bonding A(N6) with just residue information.
a_n6_single_res = (a_n6_single_h_bond.loc[:, ["eq_class_members", "don_resn", "don_resi", "don_chain"]]
                   .drop_duplicates()).rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"})
a_n6_dual_res = (a_n6_dual_h_bond.loc[:, ["eq_class_members", "don_resn", "don_resi", "don_chain"]]
                 .drop_duplicates()).rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"})
a_n6_no_res = (pd.concat([a_n6_single_res,
                         a_n6_dual_res,
                         nuc_data.loc[:, ["eq_class_members", "resn", "resi", "chain"]]])
               .drop_duplicates(keep=False))

# Prepare dataframes with just residue information of A residues that accept H-bonds at their N1, N3, or N7.
a_n1_acc_res = (acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]].eq(["N1", "A"]).all(axis='columns')]
                .loc[:, ["eq_class_members", "acc_resn", "acc_resi", "acc_chain"]].drop_duplicates()
                .rename(columns={"acc_resn": "resn", "acc_resi": "resi", "acc_chain": "chain"}))
a_n3_acc_res = (acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]].eq(["N3", "A"]).all(axis='columns')]
                .loc[:, ["eq_class_members", "acc_resn", "acc_resi", "acc_chain"]].drop_duplicates()
                .rename(columns={"acc_resn": "resn", "acc_resi": "resi", "acc_chain": "chain"}))
a_n7_acc_res = (acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]].eq(["N7", "A"]).all(axis='columns')]
                .loc[:, ["eq_class_members", "acc_resn", "acc_resi", "acc_chain"]].drop_duplicates()
                .rename(columns={"acc_resn": "resn", "acc_resi": "resi", "acc_chain": "chain"}))

# Prepare dataframes of no, single, and dual H-bonding A(N6) that also accept at the N1 with just residue information.
a_n6_no_a_n1_acc_res = a_n6_no_res.merge(a_n1_acc_res, how='inner')
a_n6_single_a_n1_acc_res = a_n6_single_res.merge(a_n1_acc_res, how='inner')
a_n6_dual_a_n1_acc_res = a_n6_dual_res.merge(a_n1_acc_res, how='inner')

# Prepare dataframes of no, single, and dual H-bonding A(N6) that also accept at the N3 with just residue information.
a_n6_no_a_n3_acc_res = a_n6_no_res.merge(a_n3_acc_res, how='inner')
a_n6_single_a_n3_acc_res = a_n6_single_res.merge(a_n3_acc_res, how='inner')
a_n6_dual_a_n3_acc_res = a_n6_dual_res.merge(a_n3_acc_res, how='inner')

# Prepare dataframes of no, single, and dual H-bonding A(N6) that also accept at the N7 with just residue information.
a_n6_no_a_n7_acc_res = a_n6_no_res.merge(a_n7_acc_res, how='inner')
a_n6_single_a_n7_acc_res = a_n6_single_res.merge(a_n7_acc_res, how='inner')
a_n6_dual_a_n7_acc_res = a_n6_dual_res.merge(a_n7_acc_res, how='inner')

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at N1 and which include H-bond information
# related to the N6 H-bond donation.
a_n6_single_a_n1_acc_h_bond_from_n6 = (a_n6_single_a_n1_acc_res
                                       .rename(columns={"resn": "don_resn", "resi": "don_resi", "chain": "don_chain"})
                                       .merge(a_n6_single_h_bond, how='inner'))
a_n6_dual_a_n1_acc_h_bond_from_n6 = (a_n6_dual_a_n1_acc_res
                                     .rename(columns={"resn": "don_resn", "resi": "don_resi", "chain": "don_chain"})
                                     .merge(a_n6_dual_h_bond, how='inner'))

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at N3 and which include H-bond information
# related to the N6 H-bond donation.
a_n6_single_a_n3_acc_h_bond_from_n6 = (a_n6_single_a_n3_acc_res
                                       .rename(columns={"resn": "don_resn", "resi": "don_resi", "chain": "don_chain"})
                                       .merge(a_n6_single_h_bond, how='inner'))
a_n6_dual_a_n3_acc_h_bond_from_n6 = (a_n6_dual_a_n3_acc_res
                                     .rename(columns={"resn": "don_resn", "resi": "don_resi", "chain": "don_chain"})
                                     .merge(a_n6_dual_h_bond, how='inner'))

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at N7 and which include H-bond information
# related to the N6 H-bond donation.
a_n6_single_a_n7_acc_h_bond_from_n6 = (a_n6_single_a_n7_acc_res
                                       .rename(columns={"resn": "don_resn", "resi": "don_resi", "chain": "don_chain"})
                                       .merge(a_n6_single_h_bond, how='inner'))
a_n6_dual_a_n7_acc_h_bond_from_n6 = (a_n6_dual_a_n7_acc_res
                                     .rename(columns={"resn": "don_resn", "resi": "don_resi", "chain": "don_chain"})
                                     .merge(a_n6_dual_h_bond, how='inner'))

# Prepare dataframes of no, single and dual A(N6) donors that also accept via the N1 including H-bond information
# related to the N1 H-bond acceptation.
a_n6_no_a_n1_acc_h_bond_to_n1 = (a_n6_no_a_n1_acc_res
                                 .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                 .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]].eq(["N1", "A"])
                                        .all(axis='columns')], how='inner'))
a_n6_single_a_n1_acc_h_bond_to_n1 = (a_n6_single_a_n1_acc_res
                                     .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                     .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]].eq(["N1", "A"])
                                            .all(axis='columns')], how='inner'))
a_n6_dual_a_n1_acc_h_bond_to_n1 = (a_n6_dual_a_n1_acc_res
                                   .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                   .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]].eq(["N1", "A"])
                                          .all(axis='columns')], how='inner'))

# Prepare dataframes of no, single and dual A(N6) donors that also accept via the N3 including H-bond information
# related to the N3 H-bond acceptation.
a_n6_no_a_n3_acc_h_bond_to_n3 = (a_n6_no_a_n3_acc_res
                                 .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                 .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]].eq(["N3", "A"])
                                        .all(axis='columns')], how='inner'))
a_n6_single_a_n3_acc_h_bond_to_n3 = (a_n6_single_a_n3_acc_res
                                     .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                     .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]].eq(["N3", "A"])
                                            .all(axis='columns')], how='inner'))
a_n6_dual_a_n3_acc_h_bond_to_n3 = (a_n6_dual_a_n3_acc_res
                                   .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                   .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]].eq(["N3", "A"])
                                          .all(axis='columns')], how='inner'))

# Prepare dataframes of no, single and dual A(N6) donors that also accept via the N7 including H-bond information
# related to the N7 H-bond acceptation.
a_n6_no_a_n7_acc_h_bond_to_n7 = (a_n6_no_a_n7_acc_res
                                 .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                 .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]].eq(["N7", "A"])
                                        .all(axis='columns')], how='inner'))
a_n6_single_a_n7_acc_h_bond_to_n7 = (a_n6_single_a_n7_acc_res
                                     .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                     .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]].eq(["N7", "A"])
                                            .all(axis='columns')], how='inner'))
a_n6_dual_a_n7_acc_h_bond_to_n7 = (a_n6_dual_a_n7_acc_res
                                   .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                   .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]].eq(["N7", "A"])
                                          .all(axis='columns')], how='inner'))

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at N1 and which include H-bond information
# related to the N6 H-bond donation and the N1 H-bond acceptation.
a_n6_single_a_n1_acc_match = (a_n6_single_a_n1_acc_h_bond_from_n6
                              .merge(a_n6_single_a_n1_acc_h_bond_to_n1,
                                     left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                     right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"], how='inner'))
a_n6_dual_a_n1_acc_match = (a_n6_dual_a_n1_acc_h_bond_from_n6
                            .merge(a_n6_dual_a_n1_acc_h_bond_to_n1,
                                   left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                   right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"], how='inner'))

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at N3 and which include H-bond information
# related to the N6 H-bond donation and the N3 H-bond acceptation.
a_n6_single_a_n3_acc_match = (a_n6_single_a_n3_acc_h_bond_from_n6
                              .merge(a_n6_single_a_n3_acc_h_bond_to_n3,
                                     left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                     right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"], how='inner'))
a_n6_dual_a_n3_acc_match = (a_n6_dual_a_n3_acc_h_bond_from_n6
                            .merge(a_n6_dual_a_n3_acc_h_bond_to_n3,
                                   left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                   right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"], how='inner'))

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at N7 and which include H-bond information
# related to the N6 H-bond donation and the N7 H-bond acceptation.
a_n6_single_a_n7_acc_match = (a_n6_single_a_n7_acc_h_bond_from_n6
                              .merge(a_n6_single_a_n7_acc_h_bond_to_n7,
                                     left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                     right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"], how='inner'))
a_n6_dual_a_n7_acc_match = (a_n6_dual_a_n7_acc_h_bond_from_n6
                            .merge(a_n6_dual_a_n7_acc_h_bond_to_n7,
                                   left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                   right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"], how='inner'))

# Prepare dataframes with just residue information of single and dual H-bonding A(N6) that also accepts an H-bond at the
# N1, N3, or N7. Additionally, there is overlap in the identity of the partner nucleobase, amino acid residue backbone,
# or amino acid residue side chain, where the N6 donates an H-bond and the N1, N3, or N7 accepts an H-bond from the same
# partner entity.
a_n6_single_a_n1_acc_ol = (a_n6_single_a_n1_acc_match[
    (((a_n6_single_a_n1_acc_match["acc_resn_x"] == a_n6_single_a_n1_acc_match["don_resn_y"]) &
      (a_n6_single_a_n1_acc_match["acc_resi_x"] == a_n6_single_a_n1_acc_match["don_resi_y"]) &
      (a_n6_single_a_n1_acc_match["acc_chain_x"] == a_n6_single_a_n1_acc_match["don_chain_y"]) &
      ~(a_n6_single_a_n1_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_single_a_n1_acc_match["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_single_a_n1_acc_match["acc_resn_x"] == a_n6_single_a_n1_acc_match["don_resn_y"]) &
      (a_n6_single_a_n1_acc_match["acc_resi_x"] == a_n6_single_a_n1_acc_match["don_resi_y"]) &
      (a_n6_single_a_n1_acc_match["acc_chain_x"] == a_n6_single_a_n1_acc_match["don_chain_y"]) &
      ~(a_n6_single_a_n1_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_single_a_n1_acc_match["don_name_y"].isin(["O2'"])) &
      (a_n6_single_a_n1_acc_match["acc_name_x"].isin(["O"])) &
      (a_n6_single_a_n1_acc_match["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())
a_n6_dual_a_n1_acc_ol = (a_n6_dual_a_n1_acc_match[
    (((a_n6_dual_a_n1_acc_match["acc_resn_x"] == a_n6_dual_a_n1_acc_match["don_resn_y"]) &
      (a_n6_dual_a_n1_acc_match["acc_resi_x"] == a_n6_dual_a_n1_acc_match["don_resi_y"]) &
      (a_n6_dual_a_n1_acc_match["acc_chain_x"] == a_n6_dual_a_n1_acc_match["don_chain_y"]) &
      ~(a_n6_dual_a_n1_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_dual_a_n1_acc_match["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_dual_a_n1_acc_match["acc_resn_x"] == a_n6_dual_a_n1_acc_match["don_resn_y"]) &
      (a_n6_dual_a_n1_acc_match["acc_resi_x"] == a_n6_dual_a_n1_acc_match["don_resi_y"]) &
      (a_n6_dual_a_n1_acc_match["acc_chain_x"] == a_n6_dual_a_n1_acc_match["don_chain_y"]) &
      ~(a_n6_dual_a_n1_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_dual_a_n1_acc_match["don_name_y"].isin(["O2'"])) &
      (a_n6_dual_a_n1_acc_match["acc_name_x"].isin(["O"])) &
      (a_n6_dual_a_n1_acc_match["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())
a_n6_single_a_n3_acc_ol = (a_n6_single_a_n3_acc_match[
    (((a_n6_single_a_n3_acc_match["acc_resn_x"] == a_n6_single_a_n3_acc_match["don_resn_y"]) &
      (a_n6_single_a_n3_acc_match["acc_resi_x"] == a_n6_single_a_n3_acc_match["don_resi_y"]) &
      (a_n6_single_a_n3_acc_match["acc_chain_x"] == a_n6_single_a_n3_acc_match["don_chain_y"]) &
      ~(a_n6_single_a_n3_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_single_a_n3_acc_match["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_single_a_n3_acc_match["acc_resn_x"] == a_n6_single_a_n3_acc_match["don_resn_y"]) &
      (a_n6_single_a_n3_acc_match["acc_resi_x"] == a_n6_single_a_n3_acc_match["don_resi_y"]) &
      (a_n6_single_a_n3_acc_match["acc_chain_x"] == a_n6_single_a_n3_acc_match["don_chain_y"]) &
      ~(a_n6_single_a_n3_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_single_a_n3_acc_match["don_name_y"].isin(["O2'"])) &
      (a_n6_single_a_n3_acc_match["acc_name_x"].isin(["O"])) &
      (a_n6_single_a_n3_acc_match["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())
a_n6_dual_a_n3_acc_ol = (a_n6_dual_a_n3_acc_match[
    (((a_n6_dual_a_n3_acc_match["acc_resn_x"] == a_n6_dual_a_n3_acc_match["don_resn_y"]) &
      (a_n6_dual_a_n3_acc_match["acc_resi_x"] == a_n6_dual_a_n3_acc_match["don_resi_y"]) &
      (a_n6_dual_a_n3_acc_match["acc_chain_x"] == a_n6_dual_a_n3_acc_match["don_chain_y"]) &
      ~(a_n6_dual_a_n3_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_dual_a_n3_acc_match["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_dual_a_n3_acc_match["acc_resn_x"] == a_n6_dual_a_n3_acc_match["don_resn_y"]) &
      (a_n6_dual_a_n3_acc_match["acc_resi_x"] == a_n6_dual_a_n3_acc_match["don_resi_y"]) &
      (a_n6_dual_a_n3_acc_match["acc_chain_x"] == a_n6_dual_a_n3_acc_match["don_chain_y"]) &
      ~(a_n6_dual_a_n3_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_dual_a_n3_acc_match["don_name_y"].isin(["O2'"])) &
      (a_n6_dual_a_n3_acc_match["acc_name_x"].isin(["O"])) &
      (a_n6_dual_a_n3_acc_match["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())
a_n6_single_a_n7_acc_ol = (a_n6_single_a_n7_acc_match[
    (((a_n6_single_a_n7_acc_match["acc_resn_x"] == a_n6_single_a_n7_acc_match["don_resn_y"]) &
      (a_n6_single_a_n7_acc_match["acc_resi_x"] == a_n6_single_a_n7_acc_match["don_resi_y"]) &
      (a_n6_single_a_n7_acc_match["acc_chain_x"] == a_n6_single_a_n7_acc_match["don_chain_y"]) &
      ~(a_n6_single_a_n7_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_single_a_n7_acc_match["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_single_a_n7_acc_match["acc_resn_x"] == a_n6_single_a_n7_acc_match["don_resn_y"]) &
      (a_n6_single_a_n7_acc_match["acc_resi_x"] == a_n6_single_a_n7_acc_match["don_resi_y"]) &
      (a_n6_single_a_n7_acc_match["acc_chain_x"] == a_n6_single_a_n7_acc_match["don_chain_y"]) &
      ~(a_n6_single_a_n7_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_single_a_n7_acc_match["don_name_y"].isin(["O2'"])) &
      (a_n6_single_a_n7_acc_match["acc_name_x"].isin(["O"])) &
      (a_n6_single_a_n7_acc_match["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())
a_n6_dual_a_n7_acc_ol = (a_n6_dual_a_n7_acc_match[
    (((a_n6_dual_a_n7_acc_match["acc_resn_x"] == a_n6_dual_a_n7_acc_match["don_resn_y"]) &
      (a_n6_dual_a_n7_acc_match["acc_resi_x"] == a_n6_dual_a_n7_acc_match["don_resi_y"]) &
      (a_n6_dual_a_n7_acc_match["acc_chain_x"] == a_n6_dual_a_n7_acc_match["don_chain_y"]) &
      ~(a_n6_dual_a_n7_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_dual_a_n7_acc_match["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_dual_a_n7_acc_match["acc_resn_x"] == a_n6_dual_a_n7_acc_match["don_resn_y"]) &
      (a_n6_dual_a_n7_acc_match["acc_resi_x"] == a_n6_dual_a_n7_acc_match["don_resi_y"]) &
      (a_n6_dual_a_n7_acc_match["acc_chain_x"] == a_n6_dual_a_n7_acc_match["don_chain_y"]) &
      ~(a_n6_dual_a_n7_acc_match["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_dual_a_n7_acc_match["don_name_y"].isin(["O2'"])) &
      (a_n6_dual_a_n7_acc_match["acc_name_x"].isin(["O"])) &
      (a_n6_dual_a_n7_acc_match["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())

# Prepare dataframes with just residue information of single and dual H-bonding A(N6) that also accepts an H-bond at the
# N1, N3, or N7. Additionally, there is no overlap in the identity of the partner nucleobase, amino acid residue
# backbone, or amino acid residue side chain. In other words, the acceptor for the N6 H-bond donation and the donor for
# the N1, N3, or N7 H-bond acceptation belong to different partner entities.
a_n6_single_a_n1_acc_no_ol = pd.concat([a_n6_single_a_n1_acc_res, a_n6_single_a_n1_acc_ol]).drop_duplicates(keep=False)
a_n6_dual_a_n1_acc_no_ol = pd.concat([a_n6_dual_a_n1_acc_res, a_n6_dual_a_n1_acc_ol]).drop_duplicates(keep=False)
a_n6_single_a_n3_acc_no_ol = pd.concat([a_n6_single_a_n3_acc_res, a_n6_single_a_n3_acc_ol]).drop_duplicates(keep=False)
a_n6_dual_a_n3_acc_no_ol = pd.concat([a_n6_dual_a_n3_acc_res, a_n6_dual_a_n3_acc_ol]).drop_duplicates(keep=False)
a_n6_single_a_n7_acc_no_ol = pd.concat([a_n6_single_a_n7_acc_res, a_n6_single_a_n7_acc_ol]).drop_duplicates(keep=False)
a_n6_dual_a_n7_acc_no_ol = pd.concat([a_n6_dual_a_n7_acc_res, a_n6_dual_a_n7_acc_ol]).drop_duplicates(keep=False)

# Prepare dataframes of A(N6) donors that also accept via the N1, N3, or N7 while only considering partner entities that
# do not overlap and which includes H-bond information related to the N1, N3, or N7 H-bond acceptation.
a_n6_single_a_n1_acc_no_ol_h_bond_to_n1 = (a_n6_single_a_n1_acc_no_ol
                                           .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                            "chain": "acc_chain"})
                                           .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]]
                                                  .eq(["N1", "A"])
                                                  .all(axis='columns')], how='inner'))
a_n6_dual_a_n1_acc_no_ol_h_bond_to_n1 = (a_n6_dual_a_n1_acc_no_ol
                                         .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                         .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]].eq(["N1", "A"])
                                                .all(axis='columns')], how='inner'))
a_n6_single_a_n3_acc_no_ol_h_bond_to_n3 = (a_n6_single_a_n3_acc_no_ol
                                           .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                            "chain": "acc_chain"})
                                           .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]]
                                                  .eq(["N3", "A"])
                                                  .all(axis='columns')], how='inner'))
a_n6_dual_a_n3_acc_no_ol_h_bond_to_n3 = (a_n6_dual_a_n3_acc_no_ol
                                         .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                         .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]].eq(["N3", "A"])
                                                .all(axis='columns')], how='inner'))
a_n6_single_a_n7_acc_no_ol_h_bond_to_n7 = (a_n6_single_a_n7_acc_no_ol
                                           .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                            "chain": "acc_chain"})
                                           .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]]
                                                  .eq(["N7", "A"])
                                                  .all(axis='columns')], how='inner'))
a_n6_dual_a_n7_acc_no_ol_h_bond_to_n7 = (a_n6_dual_a_n7_acc_no_ol
                                         .rename(columns={"resn": "acc_resn", "resi": "acc_resi", "chain": "acc_chain"})
                                         .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]].eq(["N7", "A"])
                                                .all(axis='columns')], how='inner'))

# Prepare dataframes of single and dual H-bonding A(N6) residues that donate to at least one acceptor that bears a high
# negative charge and that accept an H-bond at the N1. The acceptors could be the OD1 or OD2 of Asp, the OE1 or OE2 of
# Glu, or the OP1 or OP2 of nucleic acids. The dataframes include H-bond information related to the N6 H-bond donation.
find_neg_single_n1 = a_n6_single_a_n1_acc_h_bond_from_n6
find_neg_single_n1["neg_acc"] = pd.Series((a_n6_single_a_n1_acc_h_bond_from_n6["acc_resn"].isin(neg_acc_resn)) &
                                          (a_n6_single_a_n1_acc_h_bond_from_n6["acc_name"].isin(neg_acc_name)))
a_n6_single_a_n1_acc_h_bond_from_n6_neg = (find_neg_single_n1[find_neg_single_n1
                                           .groupby(["don_index", "eq_class_members"])["neg_acc"]
                                           .transform("any")].drop(columns=["neg_acc"]))
find_neg_dual_n1 = a_n6_dual_a_n1_acc_h_bond_from_n6
find_neg_dual_n1["neg_acc"] = pd.Series((a_n6_dual_a_n1_acc_h_bond_from_n6["acc_resn"].isin(neg_acc_resn)) &
                                        (a_n6_dual_a_n1_acc_h_bond_from_n6["acc_name"].isin(neg_acc_name)))
a_n6_dual_a_n1_acc_h_bond_from_n6_neg = (find_neg_dual_n1[find_neg_dual_n1
                                         .groupby(["don_index", "eq_class_members"])["neg_acc"].transform("any")]
                                         .drop(columns=["neg_acc"]))

# Prepare dataframes of single and dual H-bonding A(N6) residues that donate to at least one acceptor that bears a high
# negative charge and that accept an H-bond at the N3. The acceptors could be the OD1 or OD2 of Asp, the OE1 or OE2 of
# Glu, or the OP1 or OP2 of nucleic acids. The dataframes include H-bond information related to the N6 H-bond donation.
find_neg_single_n3 = a_n6_single_a_n3_acc_h_bond_from_n6
find_neg_single_n3["neg_acc"] = pd.Series((a_n6_single_a_n3_acc_h_bond_from_n6["acc_resn"].isin(neg_acc_resn)) &
                                          (a_n6_single_a_n3_acc_h_bond_from_n6["acc_name"].isin(neg_acc_name)))
a_n6_single_a_n3_acc_h_bond_from_n6_neg = (find_neg_single_n3[find_neg_single_n3
                                           .groupby(["don_index", "eq_class_members"])["neg_acc"]
                                           .transform("any")].drop(columns=["neg_acc"]))
find_neg_dual_n3 = a_n6_dual_a_n3_acc_h_bond_from_n6
find_neg_dual_n3["neg_acc"] = pd.Series((a_n6_dual_a_n3_acc_h_bond_from_n6["acc_resn"].isin(neg_acc_resn)) &
                                        (a_n6_dual_a_n3_acc_h_bond_from_n6["acc_name"].isin(neg_acc_name)))
a_n6_dual_a_n3_acc_h_bond_from_n6_neg = (find_neg_dual_n3[find_neg_dual_n3
                                         .groupby(["don_index", "eq_class_members"])["neg_acc"].transform("any")]
                                         .drop(columns=["neg_acc"]))

# Prepare dataframes of single and dual H-bonding A(N6) residues that donate to at least one acceptor that bears a high
# negative charge and that accept an H-bond at the N7. The acceptors could be the OD1 or OD2 of Asp, the OE1 or OE2 of
# Glu, or the OP1 or OP2 of nucleic acids. The dataframes include H-bond information related to the N6 H-bond donation.
find_neg_single_n7 = a_n6_single_a_n7_acc_h_bond_from_n6
find_neg_single_n7["neg_acc"] = pd.Series((a_n6_single_a_n7_acc_h_bond_from_n6["acc_resn"].isin(neg_acc_resn)) &
                                          (a_n6_single_a_n7_acc_h_bond_from_n6["acc_name"].isin(neg_acc_name)))
a_n6_single_a_n7_acc_h_bond_from_n6_neg = (find_neg_single_n7[find_neg_single_n7
                                           .groupby(["don_index", "eq_class_members"])["neg_acc"]
                                           .transform("any")].drop(columns=["neg_acc"]))
find_neg_dual_n7 = a_n6_dual_a_n7_acc_h_bond_from_n6
find_neg_dual_n7["neg_acc"] = pd.Series((a_n6_dual_a_n7_acc_h_bond_from_n6["acc_resn"].isin(neg_acc_resn)) &
                                        (a_n6_dual_a_n7_acc_h_bond_from_n6["acc_name"].isin(neg_acc_name)))
a_n6_dual_a_n7_acc_h_bond_from_n6_neg = (find_neg_dual_n7[find_neg_dual_n7
                                         .groupby(["don_index", "eq_class_members"])["neg_acc"].transform("any")]
                                         .drop(columns=["neg_acc"]))

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at the N1 with just residue information. The N6
# donates to at least one OD1 or OD2 of Asp, OE1 or OE2 of Glu, or OP1 or OP2 of nucleic acids.
a_n6_single_a_n1_acc_res_neg = ((a_n6_single_a_n1_acc_h_bond_from_n6_neg
                                .loc[:, ["eq_class_members", "don_resn", "don_resi", "don_chain"]]
                                .drop_duplicates())
                                .rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"}))
a_n6_dual_a_n1_acc_res_neg = ((a_n6_dual_a_n1_acc_h_bond_from_n6_neg
                              .loc[:, ["eq_class_members", "don_resn", "don_resi", "don_chain"]]
                              .drop_duplicates())
                              .rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"}))

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at the N3 with just residue information. The N6
# donates to at least one OD1 or OD2 of Asp, OE1 or OE2 of Glu, or OP1 or OP2 of nucleic acids.
a_n6_single_a_n3_acc_res_neg = ((a_n6_single_a_n3_acc_h_bond_from_n6_neg
                                .loc[:, ["eq_class_members", "don_resn", "don_resi", "don_chain"]]
                                .drop_duplicates())
                                .rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"}))
a_n6_dual_a_n3_acc_res_neg = ((a_n6_dual_a_n3_acc_h_bond_from_n6_neg
                              .loc[:, ["eq_class_members", "don_resn", "don_resi", "don_chain"]]
                              .drop_duplicates())
                              .rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"}))

# Prepare dataframes of single and dual H-bonding A(N6) that also accept at the N7 with just residue information. The N6
# donates to at least one OD1 or OD2 of Asp, OE1 or OE2 of Glu, or OP1 or OP2 of nucleic acids.
a_n6_single_a_n7_acc_res_neg = ((a_n6_single_a_n7_acc_h_bond_from_n6_neg
                                .loc[:, ["eq_class_members", "don_resn", "don_resi", "don_chain"]]
                                .drop_duplicates())
                                .rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"}))
a_n6_dual_a_n7_acc_res_neg = ((a_n6_dual_a_n7_acc_h_bond_from_n6_neg
                              .loc[:, ["eq_class_members", "don_resn", "don_resi", "don_chain"]]
                              .drop_duplicates())
                              .rename(columns={"don_resn": "resn", "don_resi": "resi", "don_chain": "chain"}))

# Prepare dataframes of A(N6) donors that also accept via the N1, N3, or N7, including H-bond information related to the
# N1, N3, or N7 H-bond acceptation. The N6 donates to at least one OD1 or OD2 of Asp, OE1 or OE2 of Glu, or OP1 or OP2
# of nucleic acids.
a_n6_single_a_n1_acc_h_bond_to_n1_neg = (a_n6_single_a_n1_acc_res_neg
                                         .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                          "chain": "acc_chain"})
                                         .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]]
                                                .eq(["N1", "A"])
                                                .all(axis='columns')], how='inner'))
a_n6_dual_a_n1_acc_h_bond_to_n1_neg = (a_n6_dual_a_n1_acc_res_neg
                                       .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                        "chain": "acc_chain"})
                                       .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]]
                                              .eq(["N1", "A"])
                                              .all(axis='columns')], how='inner'))
a_n6_single_a_n3_acc_h_bond_to_n3_neg = (a_n6_single_a_n3_acc_res_neg
                                         .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                          "chain": "acc_chain"})
                                         .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]]
                                                .eq(["N3", "A"])
                                                .all(axis='columns')], how='inner'))
a_n6_dual_a_n3_acc_h_bond_to_n3_neg = (a_n6_dual_a_n3_acc_res_neg
                                       .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                        "chain": "acc_chain"})
                                       .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]]
                                              .eq(["N3", "A"])
                                              .all(axis='columns')], how='inner'))
a_n6_single_a_n7_acc_h_bond_to_n7_neg = (a_n6_single_a_n7_acc_res_neg
                                         .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                          "chain": "acc_chain"})
                                         .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]]
                                                .eq(["N7", "A"])
                                                .all(axis='columns')], how='inner'))
a_n6_dual_a_n7_acc_h_bond_to_n7_neg = (a_n6_dual_a_n7_acc_res_neg
                                       .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                        "chain": "acc_chain"})
                                       .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]]
                                              .eq(["N7", "A"])
                                              .all(axis='columns')], how='inner'))

# Prepare dataframes of single and dual H-bonding A(N6) residues that only donate to OD1 or OD2 of Asp, OE1 or OE2 of
# Glu, or OP1 or OP2 of nucleic acids and that accept an H-bond at the N1. The dataframes include H-bond information
# related to the N6 H-bond donation and the N1 H-bond acceptation.
a_n6_single_a_n1_acc_match_neg = (a_n6_single_a_n1_acc_h_bond_from_n6_neg
                                  .merge(a_n6_single_a_n1_acc_h_bond_to_n1,
                                         left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                         right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"],
                                         how='inner'))
a_n6_dual_a_n1_acc_match_neg = (a_n6_dual_a_n1_acc_h_bond_from_n6_neg
                                .merge(a_n6_dual_a_n1_acc_h_bond_to_n1,
                                       left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                       right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"],
                                       how='inner'))

# Prepare dataframes of single and dual H-bonding A(N6) residues that only donate to OD1 or OD2 of Asp, OE1 or OE2 of
# Glu, or OP1 or OP2 of nucleic acids and that accept an H-bond at the N3. The dataframes include H-bond information
# related to the N6 H-bond donation and the N3 H-bond acceptation.
a_n6_single_a_n3_acc_match_neg = (a_n6_single_a_n3_acc_h_bond_from_n6_neg
                                  .merge(a_n6_single_a_n3_acc_h_bond_to_n3,
                                         left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                         right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"],
                                         how='inner'))
a_n6_dual_a_n3_acc_match_neg = (a_n6_dual_a_n3_acc_h_bond_from_n6_neg
                                .merge(a_n6_dual_a_n3_acc_h_bond_to_n3,
                                       left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                       right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"],
                                       how='inner'))

# Prepare dataframes of single and dual H-bonding A(N6) residues that only donate to OD1 or OD2 of Asp, OE1 or OE2 of
# Glu, or OP1 or OP2 of nucleic acids and that accept an H-bond at the N7. The dataframes include H-bond information
# related to the N6 H-bond donation and the N7 H-bond acceptation.
a_n6_single_a_n7_acc_match_neg = (a_n6_single_a_n7_acc_h_bond_from_n6_neg
                                  .merge(a_n6_single_a_n7_acc_h_bond_to_n7,
                                         left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                         right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"],
                                         how='inner'))
a_n6_dual_a_n7_acc_match_neg = (a_n6_dual_a_n7_acc_h_bond_from_n6_neg
                                .merge(a_n6_dual_a_n7_acc_h_bond_to_n7,
                                       left_on=["don_resn", "don_resi", "don_chain", "eq_class_members"],
                                       right_on=["acc_resn", "acc_resi", "acc_chain", "eq_class_members"],
                                       how='inner'))

# Prepare dataframes with just residue information of single and dual H-bonding A(N6) that also accepts an H-bond at the
# N1, N3, or N7. The N6 only donates to OD1 or OD2 of Asp, OE1 or OE2 of Glu, or OP1 or OP2 of nucleic acids.
# Additionally, there is overlap in the identity of the partner nucleobase, amino acid residue backbone, or amino acid
# residue side chain, where the N6 donates an H-bond and the N1, N3, or N7 accepts an H-bond from the same partner
# entity.
a_n6_single_a_n1_acc_ol_neg = (a_n6_single_a_n1_acc_match_neg[
    (((a_n6_single_a_n1_acc_match_neg["acc_resn_x"] == a_n6_single_a_n1_acc_match_neg["don_resn_y"]) &
      (a_n6_single_a_n1_acc_match_neg["acc_resi_x"] == a_n6_single_a_n1_acc_match_neg["don_resi_y"]) &
      (a_n6_single_a_n1_acc_match_neg["acc_chain_x"] == a_n6_single_a_n1_acc_match_neg["don_chain_y"]) &
      ~(a_n6_single_a_n1_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_single_a_n1_acc_match_neg["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_single_a_n1_acc_match_neg["acc_resn_x"] == a_n6_single_a_n1_acc_match_neg["don_resn_y"]) &
      (a_n6_single_a_n1_acc_match_neg["acc_resi_x"] == a_n6_single_a_n1_acc_match_neg["don_resi_y"]) &
      (a_n6_single_a_n1_acc_match_neg["acc_chain_x"] == a_n6_single_a_n1_acc_match_neg["don_chain_y"]) &
      ~(a_n6_single_a_n1_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_single_a_n1_acc_match_neg["don_name_y"].isin(["O2'"])) &
      (a_n6_single_a_n1_acc_match_neg["acc_name_x"].isin(["O"])) &
      (a_n6_single_a_n1_acc_match_neg["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())
a_n6_dual_a_n1_acc_ol_neg = (a_n6_dual_a_n1_acc_match_neg[
    (((a_n6_dual_a_n1_acc_match_neg["acc_resn_x"] == a_n6_dual_a_n1_acc_match_neg["don_resn_y"]) &
      (a_n6_dual_a_n1_acc_match_neg["acc_resi_x"] == a_n6_dual_a_n1_acc_match_neg["don_resi_y"]) &
      (a_n6_dual_a_n1_acc_match_neg["acc_chain_x"] == a_n6_dual_a_n1_acc_match_neg["don_chain_y"]) &
      ~(a_n6_dual_a_n1_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_dual_a_n1_acc_match_neg["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_dual_a_n1_acc_match_neg["acc_resn_x"] == a_n6_dual_a_n1_acc_match_neg["don_resn_y"]) &
      (a_n6_dual_a_n1_acc_match_neg["acc_resi_x"] == a_n6_dual_a_n1_acc_match_neg["don_resi_y"]) &
      (a_n6_dual_a_n1_acc_match_neg["acc_chain_x"] == a_n6_dual_a_n1_acc_match_neg["don_chain_y"]) &
      ~(a_n6_dual_a_n1_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_dual_a_n1_acc_match_neg["don_name_y"].isin(["O2'"])) &
      (a_n6_dual_a_n1_acc_match_neg["acc_name_x"].isin(["O"])) &
      (a_n6_dual_a_n1_acc_match_neg["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())
a_n6_single_a_n3_acc_ol_neg = (a_n6_single_a_n3_acc_match_neg[
    (((a_n6_single_a_n3_acc_match_neg["acc_resn_x"] == a_n6_single_a_n3_acc_match_neg["don_resn_y"]) &
      (a_n6_single_a_n3_acc_match_neg["acc_resi_x"] == a_n6_single_a_n3_acc_match_neg["don_resi_y"]) &
      (a_n6_single_a_n3_acc_match_neg["acc_chain_x"] == a_n6_single_a_n3_acc_match_neg["don_chain_y"]) &
      ~(a_n6_single_a_n3_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_single_a_n3_acc_match_neg["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_single_a_n3_acc_match_neg["acc_resn_x"] == a_n6_single_a_n3_acc_match_neg["don_resn_y"]) &
      (a_n6_single_a_n3_acc_match_neg["acc_resi_x"] == a_n6_single_a_n3_acc_match_neg["don_resi_y"]) &
      (a_n6_single_a_n3_acc_match_neg["acc_chain_x"] == a_n6_single_a_n3_acc_match_neg["don_chain_y"]) &
      ~(a_n6_single_a_n3_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_single_a_n3_acc_match_neg["don_name_y"].isin(["O2'"])) &
      (a_n6_single_a_n3_acc_match_neg["acc_name_x"].isin(["O"])) &
      (a_n6_single_a_n3_acc_match_neg["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())
a_n6_dual_a_n3_acc_ol_neg = (a_n6_dual_a_n3_acc_match_neg[
    (((a_n6_dual_a_n3_acc_match_neg["acc_resn_x"] == a_n6_dual_a_n3_acc_match_neg["don_resn_y"]) &
      (a_n6_dual_a_n3_acc_match_neg["acc_resi_x"] == a_n6_dual_a_n3_acc_match_neg["don_resi_y"]) &
      (a_n6_dual_a_n3_acc_match_neg["acc_chain_x"] == a_n6_dual_a_n3_acc_match_neg["don_chain_y"]) &
      ~(a_n6_dual_a_n3_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_dual_a_n3_acc_match_neg["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_dual_a_n3_acc_match_neg["acc_resn_x"] == a_n6_dual_a_n3_acc_match_neg["don_resn_y"]) &
      (a_n6_dual_a_n3_acc_match_neg["acc_resi_x"] == a_n6_dual_a_n3_acc_match_neg["don_resi_y"]) &
      (a_n6_dual_a_n3_acc_match_neg["acc_chain_x"] == a_n6_dual_a_n3_acc_match_neg["don_chain_y"]) &
      ~(a_n6_dual_a_n3_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_dual_a_n3_acc_match_neg["don_name_y"].isin(["O2'"])) &
      (a_n6_dual_a_n3_acc_match_neg["acc_name_x"].isin(["O"])) &
      (a_n6_dual_a_n3_acc_match_neg["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())
a_n6_single_a_n7_acc_ol_neg = (a_n6_single_a_n7_acc_match_neg[
    (((a_n6_single_a_n7_acc_match_neg["acc_resn_x"] == a_n6_single_a_n7_acc_match_neg["don_resn_y"]) &
      (a_n6_single_a_n7_acc_match_neg["acc_resi_x"] == a_n6_single_a_n7_acc_match_neg["don_resi_y"]) &
      (a_n6_single_a_n7_acc_match_neg["acc_chain_x"] == a_n6_single_a_n7_acc_match_neg["don_chain_y"]) &
      ~(a_n6_single_a_n7_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_single_a_n7_acc_match_neg["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_single_a_n7_acc_match_neg["acc_resn_x"] == a_n6_single_a_n7_acc_match_neg["don_resn_y"]) &
      (a_n6_single_a_n7_acc_match_neg["acc_resi_x"] == a_n6_single_a_n7_acc_match_neg["don_resi_y"]) &
      (a_n6_single_a_n7_acc_match_neg["acc_chain_x"] == a_n6_single_a_n7_acc_match_neg["don_chain_y"]) &
      ~(a_n6_single_a_n7_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_single_a_n7_acc_match_neg["don_name_y"].isin(["O2'"])) &
      (a_n6_single_a_n7_acc_match_neg["acc_name_x"].isin(["O"])) &
      (a_n6_single_a_n7_acc_match_neg["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())
a_n6_dual_a_n7_acc_ol_neg = (a_n6_dual_a_n7_acc_match_neg[
    (((a_n6_dual_a_n7_acc_match_neg["acc_resn_x"] == a_n6_dual_a_n7_acc_match_neg["don_resn_y"]) &
      (a_n6_dual_a_n7_acc_match_neg["acc_resi_x"] == a_n6_dual_a_n7_acc_match_neg["don_resi_y"]) &
      (a_n6_dual_a_n7_acc_match_neg["acc_chain_x"] == a_n6_dual_a_n7_acc_match_neg["don_chain_y"]) &
      ~(a_n6_dual_a_n7_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2", "O"])) &
      ~(a_n6_dual_a_n7_acc_match_neg["don_name_y"].isin(["O2'", "N"]))) |
     ((a_n6_dual_a_n7_acc_match_neg["acc_resn_x"] == a_n6_dual_a_n7_acc_match_neg["don_resn_y"]) &
      (a_n6_dual_a_n7_acc_match_neg["acc_resi_x"] == a_n6_dual_a_n7_acc_match_neg["don_resi_y"]) &
      (a_n6_dual_a_n7_acc_match_neg["acc_chain_x"] == a_n6_dual_a_n7_acc_match_neg["don_chain_y"]) &
      ~(a_n6_dual_a_n7_acc_match_neg["acc_name_x"].isin(["O2'", "O3'", "O4'", "O5'", "OP1", "OP2"])) &
      ~(a_n6_dual_a_n7_acc_match_neg["don_name_y"].isin(["O2'"])) &
      (a_n6_dual_a_n7_acc_match_neg["acc_name_x"].isin(["O"])) &
      (a_n6_dual_a_n7_acc_match_neg["don_name_y"].isin(["N"]))))]
                           .loc[:, ["eq_class_members", "don_resn_x", "don_resi_x", "don_chain_x"]]
                           .rename(columns={"don_resn_x": "resn", "don_resi_x": "resi", "don_chain_x": "chain"})
                           .drop_duplicates())

# Prepare dataframes with just residue information of single and dual H-bonding A(N6) that also accepts an H-bond at the
# N1, N3, or N7. The N6 only donates to OD1 or OD2 of Asp, OE1 or OE2 of Glu, or OP1 or OP2 of nucleic acids.
# Additionally, there is no overlap in the identity of the partner nucleobase, amino acid residue backbone, or amino
# acid residue side chain. In other words, the acceptor for the N6 H-bond donation and the donor for the N1, N3, or N7
# H-bond acceptation belong to different partner entities.
a_n6_single_a_n1_acc_no_ol_neg = (pd.concat([a_n6_single_a_n1_acc_res_neg, a_n6_single_a_n1_acc_ol_neg])
                                  .drop_duplicates(keep=False))
a_n6_dual_a_n1_acc_no_ol_neg = (pd.concat([a_n6_dual_a_n1_acc_res_neg, a_n6_dual_a_n1_acc_ol_neg])
                                .drop_duplicates(keep=False))
a_n6_single_a_n3_acc_no_ol_neg = (pd.concat([a_n6_single_a_n3_acc_res_neg, a_n6_single_a_n3_acc_ol_neg])
                                  .drop_duplicates(keep=False))
a_n6_dual_a_n3_acc_no_ol_neg = (pd.concat([a_n6_dual_a_n3_acc_res_neg, a_n6_dual_a_n3_acc_ol_neg])
                                .drop_duplicates(keep=False))
a_n6_single_a_n7_acc_no_ol_neg = (pd.concat([a_n6_single_a_n7_acc_res_neg, a_n6_single_a_n7_acc_ol_neg])
                                  .drop_duplicates(keep=False))
a_n6_dual_a_n7_acc_no_ol_neg = (pd.concat([a_n6_dual_a_n7_acc_res_neg, a_n6_dual_a_n7_acc_ol_neg])
                                .drop_duplicates(keep=False))

# Prepare dataframes of A(N6) donors that also accept via the N1, N3, or N7 while only considering partner entities that
# do not overlap and which includes H-bond information related to the N1, N3, or N7 H-bond acceptation. The N6 only
# donates to OD1 or OD2 of Asp, OE1 or OE2 of Glu, or OP1 or OP2 of nucleic acids.
a_n6_single_a_n1_acc_no_ol_h_bond_to_n1_neg = (a_n6_single_a_n1_acc_no_ol_neg
                                               .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                                "chain": "acc_chain"})
                                               .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]]
                                                      .eq(["N1", "A"])
                                                      .all(axis='columns')], how='inner'))
a_n6_dual_a_n1_acc_no_ol_h_bond_to_n1_neg = (a_n6_dual_a_n1_acc_no_ol_neg
                                             .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                              "chain": "acc_chain"})
                                             .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]]
                                                    .eq(["N1", "A"])
                                                    .all(axis='columns')], how='inner'))
a_n6_single_a_n3_acc_no_ol_h_bond_to_n3_neg = (a_n6_single_a_n3_acc_no_ol_neg
                                               .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                                "chain": "acc_chain"})
                                               .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]]
                                                      .eq(["N3", "A"])
                                                      .all(axis='columns')], how='inner'))
a_n6_dual_a_n3_acc_no_ol_h_bond_to_n3_neg = (a_n6_dual_a_n3_acc_no_ol_neg
                                             .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                              "chain": "acc_chain"})
                                             .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]]
                                                    .eq(["N3", "A"])
                                                    .all(axis='columns')], how='inner'))
a_n6_single_a_n7_acc_no_ol_h_bond_to_n7_neg = (a_n6_single_a_n7_acc_no_ol_neg
                                               .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                                "chain": "acc_chain"})
                                               .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]]
                                                      .eq(["N7", "A"])
                                                      .all(axis='columns')], how='inner'))
a_n6_dual_a_n7_acc_no_ol_h_bond_to_n7_neg = (a_n6_dual_a_n7_acc_no_ol_neg
                                             .rename(columns={"resn": "acc_resn", "resi": "acc_resi",
                                                              "chain": "acc_chain"})
                                             .merge(acc_h_bonds[acc_h_bonds[["acc_name", "acc_resn"]]
                                                    .eq(["N7", "A"])
                                                    .all(axis='columns')], how='inner'))

# Write H-bond count data regarding A(N6) donors that also accept via the N1, N3, or N7.
a_n6_donation = pd.DataFrame({
    "atom": (["A(N1)"] * 3 + ["A(N3)"] * 3 + ["A(N7)"] * 3) * 2,
    "type": [0, 1, 2] * 6,
    "overlap": ["all"] * 9 + ["no"] * 9,
    "count": [a_n6_no_a_n1_acc_res.shape[0], a_n6_single_a_n1_acc_res.shape[0], a_n6_dual_a_n1_acc_res.shape[0],
              a_n6_no_a_n3_acc_res.shape[0], a_n6_single_a_n3_acc_res.shape[0], a_n6_dual_a_n3_acc_res.shape[0],
              a_n6_no_a_n7_acc_res.shape[0], a_n6_single_a_n7_acc_res.shape[0], a_n6_dual_a_n7_acc_res.shape[0],
              a_n6_no_a_n1_acc_res.shape[0], a_n6_single_a_n1_acc_no_ol.shape[0], a_n6_dual_a_n1_acc_no_ol.shape[0],
              a_n6_no_a_n3_acc_res.shape[0], a_n6_single_a_n3_acc_no_ol.shape[0], a_n6_dual_a_n3_acc_no_ol.shape[0],
              a_n6_no_a_n7_acc_res.shape[0], a_n6_single_a_n7_acc_no_ol.shape[0], a_n6_dual_a_n7_acc_no_ol.shape[0]],
    "total": [a_n6_no_res.shape[0], a_n6_single_res.shape[0], a_n6_dual_res.shape[0]] * 6
})
a_n6_donation.to_csv(snakemake.output.h_bond_count, index=False)

# Write data on H-bond distances and angles for donors to the A(N1) that also donate via their N6 with no overlap in
# partner entity. Include data where the N6 donates to all acceptor types and where it only donates to OD1 or OD2 of
# Asp, OE1 or OE2 of Glu, or OP1 or OP2 of nucleic acids.
a_n6_don_a_n1_acc_dist_ang = pd.DataFrame({
    "overlap": (["all"] * a_n6_no_a_n1_acc_h_bond_to_n1.shape[0] +
                ["all"] * a_n6_single_a_n1_acc_h_bond_to_n1.shape[0] +
                ["all"] * a_n6_dual_a_n1_acc_h_bond_to_n1.shape[0] +
                ["all"] * a_n6_single_a_n1_acc_h_bond_to_n1_neg.shape[0] +
                ["all"] * a_n6_dual_a_n1_acc_h_bond_to_n1_neg.shape[0] +
                ["no"] * a_n6_single_a_n1_acc_no_ol_h_bond_to_n1.shape[0] +
                ["no"] * a_n6_dual_a_n1_acc_no_ol_h_bond_to_n1.shape[0] +
                ["no"] * a_n6_single_a_n1_acc_no_ol_h_bond_to_n1_neg.shape[0] +
                ["no"] * a_n6_dual_a_n1_acc_no_ol_h_bond_to_n1_neg.shape[0]),
    "type": ([0] * a_n6_no_a_n1_acc_h_bond_to_n1.shape[0] +
             [1] * a_n6_single_a_n1_acc_h_bond_to_n1.shape[0] +
             [2] * a_n6_dual_a_n1_acc_h_bond_to_n1.shape[0] +
             [1] * a_n6_single_a_n1_acc_h_bond_to_n1_neg.shape[0] +
             [2] * a_n6_dual_a_n1_acc_h_bond_to_n1_neg.shape[0] +
             [1] * a_n6_single_a_n1_acc_no_ol_h_bond_to_n1.shape[0] +
             [2] * a_n6_dual_a_n1_acc_no_ol_h_bond_to_n1.shape[0] +
             [1] * a_n6_single_a_n1_acc_no_ol_h_bond_to_n1_neg.shape[0] +
             [2] * a_n6_dual_a_n1_acc_no_ol_h_bond_to_n1_neg.shape[0]),
    "charge": (["all"] * a_n6_no_a_n1_acc_h_bond_to_n1.shape[0] +
               ["all"] * a_n6_single_a_n1_acc_h_bond_to_n1.shape[0] +
               ["all"] * a_n6_dual_a_n1_acc_h_bond_to_n1.shape[0] +
               ["neg"] * a_n6_single_a_n1_acc_h_bond_to_n1_neg.shape[0] +
               ["neg"] * a_n6_dual_a_n1_acc_h_bond_to_n1_neg.shape[0] +
               ["all"] * a_n6_single_a_n1_acc_no_ol_h_bond_to_n1.shape[0] +
               ["all"] * a_n6_dual_a_n1_acc_no_ol_h_bond_to_n1.shape[0] +
               ["neg"] * a_n6_single_a_n1_acc_no_ol_h_bond_to_n1_neg.shape[0] +
               ["neg"] * a_n6_dual_a_n1_acc_no_ol_h_bond_to_n1_neg.shape[0]),
    "h_acc_distance": (a_n6_no_a_n1_acc_h_bond_to_n1["h_acc_distance"].to_list() +
                       a_n6_single_a_n1_acc_h_bond_to_n1["h_acc_distance"].to_list() +
                       a_n6_dual_a_n1_acc_h_bond_to_n1["h_acc_distance"].to_list() +
                       a_n6_single_a_n1_acc_h_bond_to_n1_neg["h_acc_distance"].to_list() +
                       a_n6_dual_a_n1_acc_h_bond_to_n1_neg["h_acc_distance"].to_list() +
                       a_n6_single_a_n1_acc_no_ol_h_bond_to_n1["h_acc_distance"].to_list() +
                       a_n6_dual_a_n1_acc_no_ol_h_bond_to_n1["h_acc_distance"].to_list() +
                       a_n6_single_a_n1_acc_no_ol_h_bond_to_n1_neg["h_acc_distance"].to_list() +
                       a_n6_dual_a_n1_acc_no_ol_h_bond_to_n1_neg["h_acc_distance"].to_list()),
    "h_angle": (a_n6_no_a_n1_acc_h_bond_to_n1["h_angle"].to_list() +
                a_n6_single_a_n1_acc_h_bond_to_n1["h_angle"].to_list() +
                a_n6_dual_a_n1_acc_h_bond_to_n1["h_angle"].to_list() +
                a_n6_single_a_n1_acc_h_bond_to_n1_neg["h_angle"].to_list() +
                a_n6_dual_a_n1_acc_h_bond_to_n1_neg["h_angle"].to_list() +
                a_n6_single_a_n1_acc_no_ol_h_bond_to_n1["h_angle"].to_list() +
                a_n6_dual_a_n1_acc_no_ol_h_bond_to_n1["h_angle"].to_list() +
                a_n6_single_a_n1_acc_no_ol_h_bond_to_n1_neg["h_angle"].to_list() +
                a_n6_dual_a_n1_acc_no_ol_h_bond_to_n1_neg["h_angle"].to_list())
})
a_n6_don_a_n1_acc_dist_ang.to_csv(snakemake.output.a_n1_acc_dist_ang, index=False)

# Write data on H-bond distances and angles for donors to the A(N3) that also donate via their N6 with no overlap in
# partner entity. Include data where the N6 donates to all acceptor types and where it only donates to OD1 or OD2 of
# Asp, OE1 or OE2 of Glu, or OP1 or OP2 of nucleic acids.
a_n6_don_a_n3_acc_dist_ang = pd.DataFrame({
    "overlap": (["all"] * a_n6_no_a_n3_acc_h_bond_to_n3.shape[0] +
                ["all"] * a_n6_single_a_n3_acc_h_bond_to_n3.shape[0] +
                ["all"] * a_n6_dual_a_n3_acc_h_bond_to_n3.shape[0] +
                ["all"] * a_n6_single_a_n3_acc_h_bond_to_n3_neg.shape[0] +
                ["all"] * a_n6_dual_a_n3_acc_h_bond_to_n3_neg.shape[0] +
                ["no"] * a_n6_single_a_n3_acc_no_ol_h_bond_to_n3.shape[0] +
                ["no"] * a_n6_dual_a_n3_acc_no_ol_h_bond_to_n3.shape[0] +
                ["no"] * a_n6_single_a_n3_acc_no_ol_h_bond_to_n3_neg.shape[0] +
                ["no"] * a_n6_dual_a_n3_acc_no_ol_h_bond_to_n3_neg.shape[0]),
    "type": ([0] * a_n6_no_a_n3_acc_h_bond_to_n3.shape[0] +
             [1] * a_n6_single_a_n3_acc_h_bond_to_n3.shape[0] +
             [2] * a_n6_dual_a_n3_acc_h_bond_to_n3.shape[0] +
             [1] * a_n6_single_a_n3_acc_h_bond_to_n3_neg.shape[0] +
             [2] * a_n6_dual_a_n3_acc_h_bond_to_n3_neg.shape[0] +
             [1] * a_n6_single_a_n3_acc_no_ol_h_bond_to_n3.shape[0] +
             [2] * a_n6_dual_a_n3_acc_no_ol_h_bond_to_n3.shape[0] +
             [1] * a_n6_single_a_n3_acc_no_ol_h_bond_to_n3_neg.shape[0] +
             [2] * a_n6_dual_a_n3_acc_no_ol_h_bond_to_n3_neg.shape[0]),
    "charge": (["all"] * a_n6_no_a_n3_acc_h_bond_to_n3.shape[0] +
               ["all"] * a_n6_single_a_n3_acc_h_bond_to_n3.shape[0] +
               ["all"] * a_n6_dual_a_n3_acc_h_bond_to_n3.shape[0] +
               ["neg"] * a_n6_single_a_n3_acc_h_bond_to_n3_neg.shape[0] +
               ["neg"] * a_n6_dual_a_n3_acc_h_bond_to_n3_neg.shape[0] +
               ["all"] * a_n6_single_a_n3_acc_no_ol_h_bond_to_n3.shape[0] +
               ["all"] * a_n6_dual_a_n3_acc_no_ol_h_bond_to_n3.shape[0] +
               ["neg"] * a_n6_single_a_n3_acc_no_ol_h_bond_to_n3_neg.shape[0] +
               ["neg"] * a_n6_dual_a_n3_acc_no_ol_h_bond_to_n3_neg.shape[0]),
    "h_acc_distance": (a_n6_no_a_n3_acc_h_bond_to_n3["h_acc_distance"].to_list() +
                       a_n6_single_a_n3_acc_h_bond_to_n3["h_acc_distance"].to_list() +
                       a_n6_dual_a_n3_acc_h_bond_to_n3["h_acc_distance"].to_list() +
                       a_n6_single_a_n3_acc_h_bond_to_n3_neg["h_acc_distance"].to_list() +
                       a_n6_dual_a_n3_acc_h_bond_to_n3_neg["h_acc_distance"].to_list() +
                       a_n6_single_a_n3_acc_no_ol_h_bond_to_n3["h_acc_distance"].to_list() +
                       a_n6_dual_a_n3_acc_no_ol_h_bond_to_n3["h_acc_distance"].to_list() +
                       a_n6_single_a_n3_acc_no_ol_h_bond_to_n3_neg["h_acc_distance"].to_list() +
                       a_n6_dual_a_n3_acc_no_ol_h_bond_to_n3_neg["h_acc_distance"].to_list()),
    "h_angle": (a_n6_no_a_n3_acc_h_bond_to_n3["h_angle"].to_list() +
                a_n6_single_a_n3_acc_h_bond_to_n3["h_angle"].to_list() +
                a_n6_dual_a_n3_acc_h_bond_to_n3["h_angle"].to_list() +
                a_n6_single_a_n3_acc_h_bond_to_n3_neg["h_angle"].to_list() +
                a_n6_dual_a_n3_acc_h_bond_to_n3_neg["h_angle"].to_list() +
                a_n6_single_a_n3_acc_no_ol_h_bond_to_n3["h_angle"].to_list() +
                a_n6_dual_a_n3_acc_no_ol_h_bond_to_n3["h_angle"].to_list() +
                a_n6_single_a_n3_acc_no_ol_h_bond_to_n3_neg["h_angle"].to_list() +
                a_n6_dual_a_n3_acc_no_ol_h_bond_to_n3_neg["h_angle"].to_list())
})
a_n6_don_a_n3_acc_dist_ang.to_csv(snakemake.output.a_n3_acc_dist_ang, index=False)

# Write data on H-bond distances and angles for donors to the A(N7) that also donate via their N6 with no overlap in
# partner entity. Include data where the N6 donates to all acceptor types and where it only donates to OD1 or OD2 of
# Asp, OE1 or OE2 of Glu, or OP1 or OP2 of nucleic acids.
a_n6_don_a_n7_acc_dist_ang = pd.DataFrame({
    "overlap": (["all"] * a_n6_no_a_n7_acc_h_bond_to_n7.shape[0] +
                ["all"] * a_n6_single_a_n7_acc_h_bond_to_n7.shape[0] +
                ["all"] * a_n6_dual_a_n7_acc_h_bond_to_n7.shape[0] +
                ["all"] * a_n6_single_a_n7_acc_h_bond_to_n7_neg.shape[0] +
                ["all"] * a_n6_dual_a_n7_acc_h_bond_to_n7_neg.shape[0] +
                ["no"] * a_n6_single_a_n7_acc_no_ol_h_bond_to_n7.shape[0] +
                ["no"] * a_n6_dual_a_n7_acc_no_ol_h_bond_to_n7.shape[0] +
                ["no"] * a_n6_single_a_n7_acc_no_ol_h_bond_to_n7_neg.shape[0] +
                ["no"] * a_n6_dual_a_n7_acc_no_ol_h_bond_to_n7_neg.shape[0]),
    "type": ([0] * a_n6_no_a_n7_acc_h_bond_to_n7.shape[0] +
             [1] * a_n6_single_a_n7_acc_h_bond_to_n7.shape[0] +
             [2] * a_n6_dual_a_n7_acc_h_bond_to_n7.shape[0] +
             [1] * a_n6_single_a_n7_acc_h_bond_to_n7_neg.shape[0] +
             [2] * a_n6_dual_a_n7_acc_h_bond_to_n7_neg.shape[0] +
             [1] * a_n6_single_a_n7_acc_no_ol_h_bond_to_n7.shape[0] +
             [2] * a_n6_dual_a_n7_acc_no_ol_h_bond_to_n7.shape[0] +
             [1] * a_n6_single_a_n7_acc_no_ol_h_bond_to_n7_neg.shape[0] +
             [2] * a_n6_dual_a_n7_acc_no_ol_h_bond_to_n7_neg.shape[0]),
    "charge": (["all"] * a_n6_no_a_n7_acc_h_bond_to_n7.shape[0] +
               ["all"] * a_n6_single_a_n7_acc_h_bond_to_n7.shape[0] +
               ["all"] * a_n6_dual_a_n7_acc_h_bond_to_n7.shape[0] +
               ["neg"] * a_n6_single_a_n7_acc_h_bond_to_n7_neg.shape[0] +
               ["neg"] * a_n6_dual_a_n7_acc_h_bond_to_n7_neg.shape[0] +
               ["all"] * a_n6_single_a_n7_acc_no_ol_h_bond_to_n7.shape[0] +
               ["all"] * a_n6_dual_a_n7_acc_no_ol_h_bond_to_n7.shape[0] +
               ["neg"] * a_n6_single_a_n7_acc_no_ol_h_bond_to_n7_neg.shape[0] +
               ["neg"] * a_n6_dual_a_n7_acc_no_ol_h_bond_to_n7_neg.shape[0]),
    "h_acc_distance": (a_n6_no_a_n7_acc_h_bond_to_n7["h_acc_distance"].to_list() +
                       a_n6_single_a_n7_acc_h_bond_to_n7["h_acc_distance"].to_list() +
                       a_n6_dual_a_n7_acc_h_bond_to_n7["h_acc_distance"].to_list() +
                       a_n6_single_a_n7_acc_h_bond_to_n7_neg["h_acc_distance"].to_list() +
                       a_n6_dual_a_n7_acc_h_bond_to_n7_neg["h_acc_distance"].to_list() +
                       a_n6_single_a_n7_acc_no_ol_h_bond_to_n7["h_acc_distance"].to_list() +
                       a_n6_dual_a_n7_acc_no_ol_h_bond_to_n7["h_acc_distance"].to_list() +
                       a_n6_single_a_n7_acc_no_ol_h_bond_to_n7_neg["h_acc_distance"].to_list() +
                       a_n6_dual_a_n7_acc_no_ol_h_bond_to_n7_neg["h_acc_distance"].to_list()),
    "h_angle": (a_n6_no_a_n7_acc_h_bond_to_n7["h_angle"].to_list() +
                a_n6_single_a_n7_acc_h_bond_to_n7["h_angle"].to_list() +
                a_n6_dual_a_n7_acc_h_bond_to_n7["h_angle"].to_list() +
                a_n6_single_a_n7_acc_h_bond_to_n7_neg["h_angle"].to_list() +
                a_n6_dual_a_n7_acc_h_bond_to_n7_neg["h_angle"].to_list() +
                a_n6_single_a_n7_acc_no_ol_h_bond_to_n7["h_angle"].to_list() +
                a_n6_dual_a_n7_acc_no_ol_h_bond_to_n7["h_angle"].to_list() +
                a_n6_single_a_n7_acc_no_ol_h_bond_to_n7_neg["h_angle"].to_list() +
                a_n6_dual_a_n7_acc_no_ol_h_bond_to_n7_neg["h_angle"].to_list())
})
a_n6_don_a_n7_acc_dist_ang.to_csv(snakemake.output.a_n7_acc_dist_ang, index=False)
