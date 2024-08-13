"""
This script identifies candidates that may be suitable for wet-lab experiments.
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


# Define a function to determine whether a donor of interest donates to any residues not included in the allowed_acc
# list specified within the config file. If it does, exclude the donor of interest from further consideration by
# returning an empty data frame.
def check_acceptor(grp):
    for row in grp.itertuples():
        if row.h_bond == 1 and row.acc_resn not in ALLOWED_ACC:
            grp = grp.truncate(after=-1)
            return grp
    return grp


# Set the candidate criteria.
RES_MAX = snakemake.config["res_max"]
AMINE_TYPE = snakemake.config["amine_type"]
B_FACTOR = snakemake.config["candidate_b_factor_cutoff"]
ALLOWED_ACC = snakemake.config["allowed_acc"]

# Initialize the candidate dataframe.
candidate_df = pd.DataFrame(columns=['don_name', 'don_resn', 'don_resi', 'don_chain', 'b_factor', 'PDB', 'model',
                                     'eq_class_member', 'num_doi'])

# Read the data from the different equivalence class members. Filter out structures and amines based on criteria
# specified within the configuration file.
for idx in range(len(snakemake.input.data)):
    try:
        df = pd.read_csv(snakemake.input.data[idx], comment="#", keep_default_na=False, na_values="NaN",
                         dtype={"don_resi": "object", "acc_resi": "object", "don_chain": "object",
                                "acc_chain": "object"})
        # If the data frame is empty, continue the loop.
        if df.empty:
            continue
        # Find the number of unique donors of interest which will correspond to the number of A's, C's, and G's in the
        # representative RNAs.
        num_doi = len(df[df["DOI"] == 1].drop_duplicates(subset=["don_index"]))
        # If the number of unique donors of interest falls outside the specified range, continue the loop.
        if not num_doi <= RES_MAX:
            continue
        # Create a new column containing the number of unique donors of interest.
        df = df.assign(num_doi=num_doi)
        # Only keep amines that are of the specified type (e.g., no, single, or dual H-bonding) and that have an average
        # nucleobase b-factor that does not exceed the specified value.
        df = df[(df["type"] == AMINE_TYPE) & (df["b_factor"] <= B_FACTOR)]
        # Only keep amines that donate H-bonds to the specified acceptors.
        df.loc[:, ~df.columns.isin(["don_index"])] = (df.groupby('don_index', group_keys=False)
                                                      .apply(check_acceptor, include_groups=False))
        # Discard redundant entries and unnecessary columns.
        drop_columns = ['don_index', 'don_segi', 'count_1', 'count_2', 'DOI', 'don_can_NA', 'acc_index', 'acc_segi',
                        'don_acc_distance', 'h_acc_distance', 'don_angle', 'h_angle', 'h_dihedral', 'h_name', 'AOI',
                        'h_bond', 'type', 'acc_charge', 'geom', 'acc_name', 'acc_resn', 'acc_resi', 'acc_chain']
        df = df.drop_duplicates(subset=["don_index"]).drop(columns=drop_columns)
        # Concatenate the data frame containing candidates from the single equivalence class member to the data frame
        # containing candidates from all members.
        if candidate_df.empty:
            candidate_df = df.copy()
        else:
            candidate_df = pd.concat([candidate_df, df])
    # Continue the loop if one of the below errors is encountered.
    except (FileNotFoundError, pd.errors.EmptyDataError):
        continue

# Write the merged data frame to a csv file.
candidate_df.to_csv(snakemake.output.candidates, index=False, na_rep='NaN')

# Close files and reset stdout and stderr.
stdout_file.close()
stderr_file.close()
sys.stdout = stdout
sys.stderr = stderr
