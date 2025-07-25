"""
This module provides a function that reads the relevant mmCIF file to collect the specified model numbers.
"""

from gemmi import cif

# Return a list of model numbers specified in the mmCIF file.
def retrieve_model(mmCIF_dir, pdb):

    # Read the mmCIF file.
    mmCIF = cif.read_file(f"{mmCIF_dir}/{pdb.lower()}.cif")

    # Initialize the block.
    block = mmCIF.sole_block()

    # Check that the _atom_site.pdbx_PDB_model_num name exists, and retrieve its values if it does.
    if block.find_values('_atom_site.pdbx_PDB_model_num'):
        model_list = block.find_values('_atom_site.pdbx_PDB_model_num')
    else:
        return "Error"

    # Convert to a set to get unique model numbers then back to a list. Return this list.
    unique_model_list = list(set(model_list))
    unique_model_list.sort()
    return unique_model_list

# Return the latest PDB entry version specified in the mmCIF file.
def retrieve_version(mmCIF_dir, pdb):

    # Read the mmCIF file.
    mmCIF = cif.read_file(f"{mmCIF_dir}/{pdb.lower()}.cif")

    # Initialize the block.
    block = mmCIF.sole_block()

    # Return information on the latest version in the revision_history.
    if block.find('_pdbx_audit_revision_history.', ['major_revision', 'minor_revision', 'revision_date']):
        last_version = block.find('_pdbx_audit_revision_history.',
                                  ['major_revision', 'minor_revision', 'revision_date'])[-1]
        return f"{last_version[0]}.{last_version[1]}", last_version[2]
    else:
        return "Error", "Error"
