"""
This module contains various functions that work together and with PyMOL to calculate H-bonding geometry measurements.
The most notable function is evaluate. This function can be used by other modules or scripts to obtain H-bonding
geometry measurements involving specified donor and acceptor atoms. The first two arguments should describe the donor
and acceptor atoms, respectively. For each atom, a tuple containing the atom index, atom name, residue name, residue
number, and chain ID should be provided in that order. The equivalence class name and a residue library should be
provided as the third and forth arguments. The function assumes 1) that hydrogens have been added to non-rotatable donor
atoms, including both endocyclic nitrogens in histidine residues and 2) that tyrosine hydroxyl hydrogens are
approximately planar with the tyrosine side chain rings. When the potential donor or acceptor atom belongs to the side
chain of an ASN, GLN, or HIS residue, hydrogens must have been added to all the nitrogen atoms on the side chain, even
if the potential H-bond does not involve those nitrogens.
"""

from pymol import cmd
from pymol import stored
import pandas as pd


def terminal_donor(donor_atom):
    # Construct lists containing the names of the canonical protein and nucleic residues.
    protein_residues = ['ALA', 'ASP', 'ASN', 'ARG', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
                        'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    nucleic_residues = ['A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT']
    # Determine whether the donor atom is at the end of a chain and is capable of donating an H-bond.
    terminal_donating_atom = []
    if cmd.count_atoms(f'not elem H and neighbor index {donor_atom[0]}') == 1 and (
            donor_atom[1] == "N" or
            donor_atom[1] == "O5'" or
            donor_atom[1] == "O3'"):
        # An amino acid N atom at the end of the chain is capable of donating an H-bond.
        if donor_atom[1] == "N" and donor_atom[2] in protein_residues:
            terminal_donating_atom = ["N", "N", True, 0, "CA"]
        # Consider O5' and O3' atoms.
        else:
            stored.bonded_atom = ''
            cmd.iterate(f'not elem H and neighbor index {donor_atom[0]}', 'stored.bonded_atom = name')
            # Only an RNA O5' atom at the 5' end of the chain is capable of donating an H-bond.
            if (donor_atom[1] == "O5'" and stored.bonded_atom == "C5'" and donor_atom[2] in
                    nucleic_residues):
                terminal_donating_atom = ["O5'", "O", True, 0, "C5'"]
            # Only an RNA O3' atom at the 3' end of the chain is capable of donating an H-bond.
            if (donor_atom[1] == "O3'" and stored.bonded_atom == "C3'" and donor_atom[2] in
                    nucleic_residues):
                # Look for a nearby phosphorus just in case PyMOL is not detecting a bond when there really should be
                # one (check out PDB ID 1HMH as an example).
                if cmd.count_atoms(f'elem P within 1.8 of index {donor_atom[0]}') == 0:
                    terminal_donating_atom = ["O3'", "O", True, 0, "C3'"]
    return terminal_donating_atom


def evaluate(donor_atom, acceptor_atom, eq_class_mem, library):
    # Initially assume that the evaluation will complete successfully.
    successful_completion = True
    # Initialize an empty list that can be used to document why an evaluation was not completed successfully.
    notes = []
    # Construct lists containing the names of the canonical protein and nucleic residues.
    protein_residues = ['ALA', 'ASP', 'ASN', 'ARG', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
                        'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    nucleic_residues = ['A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT']
    # Determine whether the donor atom is at the end of a chain and is capable of donating an H-bond.
    terminal_donating_atom = terminal_donor(donor_atom)
    # Determine whether the donor atom is listed in the residue library.
    # If found, collect other information about that atom.
    donor_res = ''
    donor_info = []
    for residue in library:
        if donor_atom[2] == residue['res']:
            donor_res = residue['res']
            for atom in residue['don']:
                if donor_atom[1] == atom[0]:
                    donor_info = atom
                    # break the loop since the atom was found
                    break
            if (terminal_donating_atom and (donor_atom[1] == "N" and donor_atom[2] in protein_residues) or
                    ((donor_atom[1] == "O5'" or donor_atom[1] == "O3'") and donor_atom[2] in
                     nucleic_residues)):
                donor_info = terminal_donating_atom
            # break the loop since the residue was found
            break
    if not donor_res:
        successful_completion = False
        notes.append(f"Error: The donor residue {donor_atom[2]}.{donor_atom[3]}.{donor_atom[4]} of {eq_class_mem} was "
                     f"not found in the residue library when evaluating a potential H-bond.")
    if not donor_info:
        successful_completion = False
        notes.append(f"Error: The donor atom {donor_atom[1]}.{donor_atom[2]}.{donor_atom[3]}.{donor_atom[4]} of "
                     f"{eq_class_mem} was not found in the residue library when evaluating a potential H-bond.")
    # Determine whether the acceptor atom is listed in the residue library.
    # If found, collect other information about that atom.
    acceptor_res = ''
    acceptor_info = []
    for residue in library:
        if acceptor_atom[2] == residue['res']:
            acceptor_res = residue['res']
            for atom in residue['acc']:
                if acceptor_atom[1] == atom[0]:
                    acceptor_info = atom
                    # break the loop since the atom was found
                    break
            # break the loop since the residue was found
            break
    if not acceptor_res:
        successful_completion = False
        notes.append(f"Error: The acceptor residue {acceptor_atom[2]}.{acceptor_atom[3]}.{acceptor_atom[4]} of "
                     f"{eq_class_mem} was not found in the residue library when evaluating a potential H-bond.")
    if not acceptor_info:
        successful_completion = False
        notes.append(f"Error: The acceptor atom {acceptor_atom[1]}.{acceptor_atom[2]}.{acceptor_atom[3]}."
                     f"{acceptor_atom[4]} of {eq_class_mem} was not found in the residue library when evaluating a "
                     f"potential H-bond.")
    # If any of the residues of the donor or acceptor atoms or the atoms themselves were not found in the residue
    # library, return a non-successful completion with the relevant notes.
    if not successful_completion:
        return [successful_completion, notes]
    # Get the distance and angle values for the donor/acceptor pair.
    geometry = calc_geom(donor_atom, acceptor_atom, donor_info, acceptor_info, eq_class_mem)
    # If the geometry calculation was not successful, return an explanation of what went wrong.
    geometry_successful = geometry[0]
    if not geometry_successful:
        return geometry
    # Proceed if the geometry calculation was successful.
    else:
        # Return the H-bonding geometry measurements.
        if len(geometry[1]) == 1:
            first_set = geometry[1][0]
            return [geometry_successful, [first_set]]
        if len(geometry[1]) == 2:
            first_set = geometry[1][0]
            second_set = geometry[1][1]
            return [geometry_successful, [first_set, second_set]]


def calc_geom(donor_atom, acceptor_atom, donor_info, acceptor_info, eq_class_mem):
    # Initially assume that the evaluation will complete successfully.
    successful_completion = True
    # Initialize an empty list that can be used to document why an evaluation was not completed successfully.
    notes = []
    # Get the D-A distance.
    don_acc_distance = cmd.get_distance(f'index {donor_atom[0]}', f'index {acceptor_atom[0]}')
    # If the donor is non-rotatable, use the donor hydrogen(s) locations to calculate the H-A distance(s) and angle(s).
    # Additionally, if the hydrogens belong to an RNA exocyclic amine, calculate the dihedral between the hydrogen
    # closest to the WCF nucleobase edge and the nearby endocyclic nitrogen that is part of the 6-membered ring. For
    # instance, this would be the N1-C6-N6-H dihedral for A residues.
    if not donor_info[2]:
        # Collect the indices of the hydrogens.
        stored.hydrogen = []
        cmd.iterate(f'elem H and neighbor index {donor_atom[0]}', 'stored.hydrogen.append((index, '
                                                                  'name, resn, resi, chain))')
        # Issue an error message if there are no hydrogens.
        if len(stored.hydrogen) == 0:
            successful_completion = False
            notes.append(f"Error: The non-rotatable donor atom {donor_atom[1]}.{donor_atom[2]}.{donor_atom[3]}."
                         f"{donor_atom[4]} of {eq_class_mem} has zero hydrogens.")
            return [successful_completion, notes]
        # Get the measurements for donors that have one hydrogen.
        elif len(stored.hydrogen) == 1:
            if donor_atom[1:3] in [["N6", "A"], ["N4", "C"], ["N2", "G"]]:
                successful_completion = False
                notes.append(f"Error: The exocyclic amine {donor_atom[1]}.{donor_atom[2]}.{donor_atom[3]}."
                             f"{donor_atom[4]} of {eq_class_mem} has only one hydrogen.")
                return [successful_completion, notes]
            else:
                h_acc_distance = [cmd.get_distance(f'index {stored.hydrogen[0][0]}', f'index {acceptor_atom[0]}'),
                                  pd.NA]
                h_angle = [cmd.get_angle(f'index {donor_atom[0]}', f'index {stored.hydrogen[0][0]}',
                                         f'index {acceptor_atom[0]}'), pd.NA]
                h_dihedral = [pd.NA, pd.NA]
                h_name = [stored.hydrogen[0][1], pd.NA]
                # Return the geometry values.
                return [successful_completion,
                        [[don_acc_distance, h_acc_distance[0], h_angle[0], h_dihedral[0],
                          h_name[0]]]]
        # Get the measurements for donors that have two hydrogens.
        elif len(stored.hydrogen) == 2:
            # Collect information about the nearby endocyclic nitrogen.
            if donor_atom[2] in ["A", "G", "DA", "DG"]:
                stored.end_n = []
                cmd.iterate(f'name N1 and neighbor index {stored.don_ant[0][0]}',
                            'stored.end_n.append((index, name, resn, resi, chain))')
            elif donor_atom[2] in ["C", "DC"]:
                stored.end_n = []
                cmd.iterate(f'name N3 and neighbor index {stored.don_ant[0][0]}',
                            'stored.end_n.append((index, name, resn, resi, chain))')
            else:
                stored.end_n = [pd.NA]
            # Issue an error message if the number of identified endocyclic nitrogen atoms does not equal one.
            if len(stored.end_n) != 1:
                successful_completion = False
                notes.append(f"Error: The number of identified endocyclic nitrogen atoms does not equal one for "
                             f"residue {donor_atom[2]}.{donor_atom[3]}.{donor_atom[4]} of {eq_class_mem}.")
                return [successful_completion, notes]
            h_acc_distance = [cmd.get_distance(f'index {stored.hydrogen[0][0]}', f'index {acceptor_atom[0]}'),
                              cmd.get_distance(f'index {stored.hydrogen[1][0]}', f'index {acceptor_atom[0]}')]
            h_angle = [cmd.get_angle(f'index {donor_atom[0]}', f'index {stored.hydrogen[0][0]}',
                                     f'index {acceptor_atom[0]}'),
                       cmd.get_angle(f'index {donor_atom[0]}', f'index {stored.hydrogen[1][0]}',
                                     f'index {acceptor_atom[0]}')]
            if not pd.isna(stored.end_n[0]):
                h_dihedral = [cmd.get_dihedral(f'index {stored.end_n[0][0]}', f'index {stored.don_ant[0][0]}',
                                               f'index {donor_atom[0]}', f'index {stored.hydrogen[0][0]}'),
                              cmd.get_dihedral(f'index {stored.end_n[0][0]}', f'index {stored.don_ant[0][0]}',
                                               f'index {donor_atom[0]}', f'index {stored.hydrogen[1][0]}')]
            else:
                h_dihedral = [pd.NA, pd.NA]
            h_name = [stored.hydrogen[0][1], stored.hydrogen[1][1]]
            # Return the geometry values.
            return [successful_completion,
                    [[don_acc_distance, h_acc_distance[0], h_angle[0], h_dihedral[0], h_name[0]],
                     [don_acc_distance, h_acc_distance[1], h_angle[1], h_dihedral[1], h_name[1]]]]
        # Issue an error message if there are more than two hydrogens.
        elif len(stored.hydrogen) > 2:
            successful_completion = False
            notes.append(f"Error: The non-rotatable donor atom {donor_atom[1]}.{donor_atom[2]}.{donor_atom[3]}."
                         f"{donor_atom[4]} of {eq_class_mem} has more than two hydrogens.")
            return [successful_completion, notes]
    else:
        h_acc_distance = [pd.NA, pd.NA]
        h_angle = [pd.NA, pd.NA]
        h_dihedral = [pd.NA, pd.NA]
        h_name = [pd.NA, pd.NA]
        # Return the geometry values.
        return [successful_completion,
                [[don_acc_distance, h_acc_distance[0], h_angle[0], h_dihedral[0], h_name[0]]]]
