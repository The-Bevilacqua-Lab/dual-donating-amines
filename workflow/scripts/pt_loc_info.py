"""
This script reads csv files corresponding to locations 1 and 2 on the dual-donating amine pseudotorsion plots. The csv
files contain information on some of the residues contributing to the counts in those locations. It then collects
information on the identity of the nearby residues. This can be used to create fragments that are extracted from the
structures which can then be analyzed using clustering.
"""

import sys
import subprocess
import pandas as pd
from Bio.PDB.MMCIFParser import MMCIFParser


# Redirect stdout and stderr to log files.
stdout = sys.stdout
stderr = sys.stderr
stdout_file = open(snakemake.log.stdout, mode='w')
stderr_file = open(snakemake.log.stderr, mode='w')
sys.stdout = stdout_file
sys.stderr = stderr_file

# The following block gathers residue information for location 1.
location_1_df = pd.read_csv(snakemake.input["location_1_residues"], keep_default_na=False, na_values="NaN", dtype="str")
location_1_frag_info = []
for row in location_1_df.itertuples():

    # Create a Biopython structure object.
    parser = MMCIFParser(QUIET=True)
    pdb_file = f'{snakemake.config["original_mmcif_dir"]}{row.PDB.lower()}.cif'
    structure = parser.get_structure(row.PDB, pdb_file)

    # Set variables pertaining to the dual donating A(N6).
    # Residue numbers and insertion codes are reported together in row.don_resi.
    # If an insertion code is present, tease it apart from the residue number.
    if row.don_resi[-1].isdigit():
        a_don_resi_int = int(row.don_resi)
        a_don_resi_icode = " "
    else:
        a_don_resi_int = int(row.don_resi[:-1])
        a_don_resi_icode = row.don_resi[-1]
        print(f"Note: One or more residues in the fragment for PDB ID {row.PDB} contains an insertion code.")
        # This script does not currently handle insertion codes that are more than one character in length.
        if not row.don_resi[-2].isdigit():
            print(f"Error: PDB ID {row.PDB} contains an insertion code that is more than one character in length.")
            continue
    a_don_chain_str = row.don_chain
    a_don_chain = structure[int(row.model) - 1][a_don_chain_str]
    a_don_chain_residues = list(a_don_chain.get_residues())
    a_don_neighbor_list = []
    for idx, residue in enumerate(a_don_chain_residues):
        if residue.get_id()[1] == a_don_resi_int and residue.get_id()[2] == a_don_resi_icode:
            a_don_neighbor_list = a_don_chain_residues[idx - 2:idx + 4]
            break

    # Set variables pertaining to the acceptors of the dual donating A(N6).
    if row.acc_pair_1_name == "N7":
        # Residue numbers and insertion codes are reported together in row.acc_pair_1_resi.
        # If an insertion code is present, tease it apart from the residue number.
        if row.acc_pair_1_resi[-1].isdigit():
            n7_resi_int = int(row.acc_pair_1_resi)
            n7_resi_icode = " "
            n7_chain_str = row.acc_pair_1_chain
        else:
            n7_resi_int = int(row.acc_pair_1_resi[:-1])
            n7_resi_icode = row.acc_pair_1_resi[-1]
            n7_chain_str = row.acc_pair_1_chain
            print(f"Note: One or more residues in the fragment for PDB ID {row.PDB} contains an insertion code.")
            # This script does not currently handle insertion codes that are more than one character in length.
            if not row.acc_pair_1_resi[-2].isdigit():
                print(f"Error: PDB ID {row.PDB} contains an insertion code that is more than one character in length.")
                continue
    else:
        # Residue numbers and insertion codes are reported together in row.acc_pair_2_resi.
        # If an insertion code is present, tease it apart from the residue number.
        if row.acc_pair_2_resi[-1].isdigit():
            n7_resi_int = int(row.acc_pair_2_resi)
            n7_resi_icode = " "
            n7_chain_str = row.acc_pair_2_chain
        else:
            n7_resi_int = int(row.acc_pair_2_resi[:-1])
            n7_resi_icode = row.acc_pair_2_resi[-1]
            n7_chain_str = row.acc_pair_2_chain
            print(f"Note: One or more residues in the fragment for PDB ID {row.PDB} contains an insertion code.")
            # This script does not currently handle insertion codes that are more than one character in length.
            if not row.acc_pair_2_resi[-2].isdigit():
                print(f"Error: PDB ID {row.PDB} contains an insertion code that is more than one character in length.")
                continue
    n7_chain = structure[int(row.model) - 1][n7_chain_str]
    n7_chain_residues = list(n7_chain.get_residues())
    n7_neighbor_list = []
    for idx, residue in enumerate(n7_chain_residues):
        if residue.get_id()[1] == n7_resi_int and residue.get_id()[2] == n7_resi_icode:
            n7_neighbor_list = n7_chain_residues[idx - 2:idx + 3]
            break

    # Append the information to the frag info list.
    a_don_neighbor_resi_list = []
    a_don_neighbor_icode_list = []
    for residue in a_don_neighbor_list:
        a_don_neighbor_resi_list.append(residue.get_id()[1])
        a_don_neighbor_icode_list.append(residue.get_id()[2])
    n7_neighbor_resi_list = []
    n7_neighbor_icode_list = []
    for residue in n7_neighbor_list:
        n7_neighbor_resi_list.append(residue.get_id()[1])
        n7_neighbor_icode_list.append(residue.get_id()[2])
    location_1_frag_info.append([row.PDB, row.model] +
                                [a_don_chain_str] + a_don_neighbor_resi_list + a_don_neighbor_icode_list +
                                [n7_chain_str] + n7_neighbor_resi_list + n7_neighbor_icode_list)

# Prepare a data frame containing the fragment info and then save it to a csv file.
location_1_frag_info_df = pd.DataFrame(location_1_frag_info, columns =
                                       ["pdb", "model",
                                       "chain_A",
                                       "resi_A_1", "resi_A_2", "resi_A_3", "resi_A_4", "resi_A_5", "resi_A_6",
                                       "icode_A_1", "icode_A_2", "icode_A_3", "icode_A_4", "icode_A_5", "icode_A_6",
                                       "chain_B",
                                       "resi_B_1", "resi_B_2", "resi_B_3", "resi_B_4", "resi_B_5",
                                       "icode_B_1", "icode_B_2", "icode_B_3", "icode_B_4", "icode_B_5"])
location_1_frag_info_df.to_csv(snakemake.output["location_1_frag_info"], index=False, na_rep='NaN')

# The following block gathers residue information for location 2.
location_2_df = pd.read_csv(snakemake.input["location_2_residues"], keep_default_na=False, na_values="NaN", dtype="str")
location_2_frag_info = []
for row in location_2_df.itertuples():

    # Create a Biopython structure object.
    parser = MMCIFParser(QUIET=True)
    pdb_file = f'{snakemake.config["original_mmcif_dir"]}{row.PDB.lower()}.cif'
    structure = parser.get_structure(row.PDB, pdb_file)

    # Set variables pertaining to the dual donating G(N2).
    # Residue numbers and insertion codes are reported together in row.don_resi.
    # If an insertion code is present, tease it apart from the residue number.
    if row.don_resi[-1].isdigit():
        g_resi_int = int(row.don_resi)
        g_resi_icode = " "
    else:
        g_resi_int = int(row.don_resi[:-1])
        g_resi_icode = row.don_resi[-1]
        print(f"Note: One or more residues in the fragment for PDB ID {row.PDB} contains an insertion code.")
        # This script does not currently handle insertion codes that are more than one character in length.
        if not row.don_resi[-2].isdigit():
            print(f"Error: PDB ID {row.PDB} contains an insertion code that is more than one character in length.")
            continue
    g_chain_str = row.don_chain
    g_chain = structure[int(row.model)-1][g_chain_str]
    g_chain_residues = list(g_chain.get_residues())
    g_neighbor_list = []
    for idx, residue in enumerate(g_chain_residues):
        if residue.get_id()[1] == g_resi_int and residue.get_id()[2] == g_resi_icode:
            g_neighbor_list = g_chain_residues[idx-3:idx+3]
            break

    # Set variables pertaining to the acceptors of the dual donating G(N2).
    if not row.acc_pair_1_name == "O4":
        # Residue numbers and insertion codes are reported together in row.acc_pair_1_resi.
        # If an insertion code is present, tease it apart from the residue number.
        if row.acc_pair_1_resi[-1].isdigit():
            npo_o5p_resi_int = int(row.acc_pair_1_resi)
            npo_o5p_resi_icode = " "
            npo_o5p_chain_str = row.acc_pair_1_chain
        else:
            npo_o5p_resi_int = int(row.acc_pair_1_resi[:-1])
            npo_o5p_resi_icode = row.acc_pair_1_resi[-1]
            npo_o5p_chain_str = row.acc_pair_1_chain
            print(f"Note: One or more residues in the fragment for PDB ID {row.PDB} contains an insertion code.")
            # This script does not currently handle insertion codes that are more than one character in length.
            if not row.acc_pair_1_resi[-2].isdigit():
                print(f"Error: PDB ID {row.PDB} contains an insertion code that is more than one character in length.")
                continue
    else:
        # Residue numbers and insertion codes are reported together in row.acc_pair_2_resi.
        # If an insertion code is present, tease it apart from the residue number.
        if row.acc_pair_2_resi[-1].isdigit():
            npo_o5p_resi_int = int(row.acc_pair_2_resi)
            npo_o5p_resi_icode = " "
            npo_o5p_chain_str = row.acc_pair_2_chain
        else:
            npo_o5p_resi_int = int(row.acc_pair_2_resi[:-1])
            npo_o5p_resi_icode = row.acc_pair_2_resi[-1]
            npo_o5p_chain_str = row.acc_pair_2_chain
            print(f"Note: One or more residues in the fragment for PDB ID {row.PDB} contains an insertion code.")
            # This script does not currently handle insertion codes that are more than one character in length.
            if not row.acc_pair_2_resi[-2].isdigit():
                print(f"Error: PDB ID {row.PDB} contains an insertion code that is more than one character in length.")
                continue
    npo_o5p_chain = structure[int(row.model)-1][npo_o5p_chain_str]
    npo_o5p_chain_residues = list(npo_o5p_chain.get_residues())
    npo_o5p_neighbor_list = []
    for idx, residue in enumerate(npo_o5p_chain_residues):
        if residue.get_id()[1] == npo_o5p_resi_int and residue.get_id()[2] == npo_o5p_resi_icode:
            npo_o5p_neighbor_list = npo_o5p_chain_residues[idx-1:idx+4]
            break

    # Append the information to the frag info list.
    g_neighbor_resi_list = []
    g_neighbor_icode_list = []
    for residue in g_neighbor_list:
        g_neighbor_resi_list.append(residue.get_id()[1])
        g_neighbor_icode_list.append(residue.get_id()[2])
    npo_o5p_neighbor_resi_list = []
    npo_o5p_neighbor_icode_list = []
    for residue in npo_o5p_neighbor_list:
        npo_o5p_neighbor_resi_list.append(residue.get_id()[1])
        npo_o5p_neighbor_icode_list.append(residue.get_id()[2])
    location_2_frag_info.append([row.PDB, row.model] +
                                [g_chain_str] + g_neighbor_resi_list + g_neighbor_icode_list +
                                [npo_o5p_chain_str] + npo_o5p_neighbor_resi_list + npo_o5p_neighbor_icode_list)

# Prepare a data frame containing the fragment info and then save it to a csv file.
location_2_frag_info_df = pd.DataFrame(location_2_frag_info, columns =
                                       ["pdb", "model",
                                        "chain_A",
                                        "resi_A_1", "resi_A_2", "resi_A_3", "resi_A_4", "resi_A_5", "resi_A_6",
                                        "icode_A_1", "icode_A_2", "icode_A_3", "icode_A_4", "icode_A_5", "icode_A_6",
                                        "chain_B",
                                        "resi_B_1", "resi_B_2", "resi_B_3", "resi_B_4", "resi_B_5",
                                        "icode_B_1", "icode_B_2", "icode_B_3", "icode_B_4", "icode_B_5"])
location_2_frag_info_df.to_csv(snakemake.output["location_2_frag_info"], index=False, na_rep='NaN')

# Remove the original mmCIF folder if specified.
if snakemake.config["remove_original_mmcif"]:
    subprocess.run(["rm", "-r", snakemake.config["original_mmcif_dir"]])

# Remove the disordered mmCIF folder if specified.
if snakemake.config["remove_disordered_mmcif"]:
    subprocess.run(["rm", "-r", snakemake.config["disordered_mmcif_dir"]])

# Close files and reset stdout and stderr.
stdout_file.close()
stderr_file.close()
sys.stdout = stdout
sys.stderr = stderr
