"""
This module reads the relevant mmCIF file to collect the specified model numbers.
"""

from gemmi import cif

# Return a list of model numbers specified in the mmCIF file.
def retrieve_model(mmCIF_dir, pdb):

    # Read the mmCIF file.
    mmCIF = cif.read_file(f"{mmCIF_dir}/{pdb.lower()}.cif")

    # Initialize the block.
    block = mmCIF.sole_block()

    # Check that the _atom_site.pdbx_PDB_model_num name exists, and retrieve its values if it does.
    model_list = []
    if block.find_values('_atom_site.pdbx_PDB_model_num'):
        model_list = block.find_values('_atom_site.pdbx_PDB_model_num')
    else:
        return "Error"

    # Convert to a set to get unique model numbers then back to a list. Return this list.
    unique_model_list = list(set(model_list))
    unique_model_list.sort()
    return unique_model_list
