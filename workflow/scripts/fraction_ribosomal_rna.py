import sys
import pandas as pd

# Redirect stdout and stderr to log files.
stdout = sys.stdout
stderr = sys.stderr
stdout_file = open(snakemake.log.stdout, mode='w')
stderr_file = open(snakemake.log.stderr, mode='w')
sys.stdout = stdout_file
sys.stderr = stderr_file


# This function is to be run if an error is encountered.
def error(msg):
    # Print the error message.
    print(msg)
    # Close files, reset stdout and stderr, and exit.
    stdout_file.close()
    stderr_file.close()
    sys.stdout = stdout
    sys.stderr = stderr
    sys.exit(0)


# This function finds the standardized name that is applicable to the chain of the amine.
def find_name(row):
    index = -1
    for idx, ife_element in enumerate(row["ife_id"].split("+")):
        chain = ife_element.split("|")[2]
        if chain == row["don_chain"]:
            index = idx
    if index == -1:
        error('Error: The standardized name of the chain could not be found for the amine with index '
              f'{row["don_index"]} from equivalence class member {row["eq_class_member"]}.')
    std_name = row["standardized_name"].split("+")[index] if not pd.isna(row["standardized_name"]) else pd.NA
    return std_name


# Read the combined data.
combined_df = pd.read_csv(snakemake.input.combined, comment="#", keep_default_na=False, na_values="NaN", dtype="string")
ifes_df = pd.read_csv(f"{snakemake.config['rep_set_file'][:-4]}_filtered.csv",
                      keep_default_na=False, na_values="NA", dtype="string")

# Create a data frame with unique rows of amines and information from the filtered IFEs file.
amines_df = combined_df.drop_duplicates(subset=["don_index", "eq_class_member"])
amines_df = amines_df.assign(ec_id=amines_df["eq_class_member"].replace(regex="([^_]*_[^_]*_[^_]*).*", value="\\1"))
amines_df = amines_df.merge(ifes_df, on="ec_id")
amines_df = amines_df.assign(don_std_name=amines_df.apply(find_name, axis=1, result_type='expand'))

# Write to the output file.
with open(snakemake.output.note, "w") as f:

    # How many IFEs contain one or more ribosomal RNAs?
    ribosomal_ifes = len(ifes_df.loc[ifes_df["standardized_name"].str.contains("ribosomal RNA")
                                     & ~ifes_df["standardized_name"].isna(), :])
    other_ifes = len(ifes_df.loc[~(ifes_df["standardized_name"].str.contains("ribosomal RNA")
                                   & ~ifes_df["standardized_name"].isna()), :])
    f.write(f"There are {ribosomal_ifes} IFEs that contain one or more ribosomal RNAs and {other_ifes} that contain no "
            "ribosomal RNAs.\n")

    # How many amines are from chains that are ribosomal RNAs?
    ribosomal_amines = len(amines_df.loc[amines_df["don_std_name"].str.contains("ribosomal RNA")
                                         & ~amines_df["don_std_name"].isna(), :])
    other_amines = len(amines_df.loc[~(amines_df["don_std_name"].str.contains("ribosomal RNA")
                                       & ~amines_df["don_std_name"].isna()), :])
    f.write(f"There are {ribosomal_amines} amines that belong to a ribosomal RNA and {other_amines} amines that belong "
            "to other RNAs.")

# Close files and reset stdout and stderr.
stdout_file.close()
stderr_file.close()
sys.stdout = stdout
sys.stderr = stderr
