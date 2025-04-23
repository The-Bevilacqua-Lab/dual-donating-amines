"""
This module contains the evaluate function which works with PyMOL to calculate H-bonding geometry measurements. This
function can be used by other modules or scripts to obtain H-bonding geometry measurements involving specified donor
and acceptor atoms, where the donor is a nucleobase amine. The first two arguments should describe the donor and
acceptor atoms, respectively. For each atom, a tuple containing the atom index, atom name, residue name, residue
number, and chain ID should be provided in that order. The equivalence class member ID and a residue library should be
provided as the third and forth arguments. The amine donor atom must be bonded to exactly two hydrogens.
"""

from pymol import cmd
from pymol import stored
import pandas as pd


def evaluate(donor_atom, acceptor_atom, eq_class_mem, library, dual_h_donors):
    # Return NaN values if the donor is not included in the config file's dual_h_donors list.
    if f"{donor_atom[2]}.{donor_atom[1]}" not in dual_h_donors:
        return [True, [[pd.NA, pd.NA, pd.NA, pd.NA, pd.NA, pd.NA]]]
    # Construct lists containing the names of the canonical protein and nucleic residues.
    protein_residues = ['ALA', 'ASP', 'ASN', 'ARG', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
                        'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    nucleic_residues = ['A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT']
    # Determine whether the donor atom is listed in the residue library.
    # If found, collect other information about that atom.
    donor_info = []
    for residue in library:
        if donor_atom[2] == residue['res']:
            for atom in residue['don']:
                if donor_atom[1] == atom[0]:
                    donor_info = atom
                    # Break the loop since the atom was found.
                    break
            # Break the loop since the residue was found.
            break
    # If any of the donor atoms were not found in the residue library, return a non-successful completion.
    if not donor_info:
        return [False, f"Error: The donor atom {donor_atom[4]}.{donor_atom[2]}.{donor_atom[3]}.{donor_atom[1]} of "
                       f"{eq_class_mem} was not found in the residue library when evaluating a potential H-bond."]
    # Get the D-A distance.
    don_acc_distance = cmd.get_distance(f'index {donor_atom[0]}', f'index {acceptor_atom[0]}')
    # Collect information about the donor antecedent atom.
    stored.don_ant = []
    cmd.iterate(f'name {donor_info[4]} and neighbor index {donor_atom[0]}',
                'stored.don_ant.append((index, name, resn, resi, chain))')
    # Issue an error message if the number of identified donor antecedent atoms does not equal one.
    if len(stored.don_ant) != 1:
        return [False, f"Error: The number of identified donor antecedent atoms does not equal one for {donor_atom[4]}."
                       f"{donor_atom[2]}.{donor_atom[3]}.{donor_atom[1]} of {eq_class_mem}."]
    # Calculate the H-A distances, D-H-A angles, and the dihedrals between the hydrogens and the WCF endocyclic
    # nitrogen. For instance, an adenine dihedral would be the N1-C6-N6-H dihedral for one of the amine hydrogens.
    # Collect the indices of the hydrogens.
    stored.hydrogen = []
    cmd.iterate(f'elem H and neighbor index {donor_atom[0]}',
                'stored.hydrogen.append((index, name, resn, resi, chain))')
    # Issue an error message if there is not exactly two hydrogens.
    if len(stored.hydrogen) != 2:
        return [False, f"Error: The amine {donor_atom[4]}.{donor_atom[2]}.{donor_atom[3]}.{donor_atom[1]} of "
                       f"{eq_class_mem} does not have exactly two hydrogens."]
    # Get the measurements for donors that have two hydrogens.
    elif f"{donor_atom[2]}.{donor_atom[1]}" in dual_h_donors:
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
            stored.end_n = []
        # Issue an error message if the number of identified endocyclic nitrogen atoms does not equal one.
        if len(stored.end_n) != 1:
            return [False, f"Error: The number of identified endocyclic nitrogen atoms does not equal one for "
                           f"residue {donor_atom[4]}.{donor_atom[2]}.{donor_atom[3]} of {eq_class_mem}."]
        h_acc_distance = [cmd.get_distance(f'index {stored.hydrogen[0][0]}', f'index {acceptor_atom[0]}'),
                          cmd.get_distance(f'index {stored.hydrogen[1][0]}', f'index {acceptor_atom[0]}')]
        h_angle = [cmd.get_angle(f'index {donor_atom[0]}', f'index {stored.hydrogen[0][0]}',
                                 f'index {acceptor_atom[0]}'),
                   cmd.get_angle(f'index {donor_atom[0]}', f'index {stored.hydrogen[1][0]}',
                                 f'index {acceptor_atom[0]}')]
        h_dihedral = [cmd.get_dihedral(f'index {stored.end_n[0][0]}', f'index {stored.don_ant[0][0]}',
                                       f'index {donor_atom[0]}', f'index {stored.hydrogen[0][0]}'),
                      cmd.get_dihedral(f'index {stored.end_n[0][0]}', f'index {stored.don_ant[0][0]}',
                                       f'index {donor_atom[0]}', f'index {stored.hydrogen[1][0]}')]
        h_name = [stored.hydrogen[0][1], stored.hydrogen[1][1]]
        # Return the geometry values.
        return [True, [[don_acc_distance, h_acc_distance[0], h_angle[0], h_dihedral[0], h_name[0]],
                       [don_acc_distance, h_acc_distance[1], h_angle[1], h_dihedral[1], h_name[1]]]]
    else:
        return [False, f"Error: The atom {donor_atom[4]}.{donor_atom[2]}.{donor_atom[3]}.{donor_atom[1]} of "
                       f"{eq_class_mem} is not included within the config file's dual_h_donors."]
