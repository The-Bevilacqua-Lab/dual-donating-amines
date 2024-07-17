"""
The remove function in this module removes alternate conformations such that each atom in the resulting structure only
has a single conformation. Alternate conformations that have the highest occupancy factor are kept while the others are
removed. If a set of alternate conformations have the same occupancy factor, the conformation with the lowest alt ID is
retained. For instance, if two conformations both have an occupancy factor of 0.5 and have alt IDs A and B, the
conformation with alt ID A would be retained. When calling this function, only one object should be loaded into the
PyMOL session.
"""

from pymol import cmd
from pymol import stored


def remove(eq_class_mem):
    # store atoms that have alternate conformations
    stored.alt_atoms = []
    cmd.iterate('not alt ""', 'stored.alt_atoms.append((name,resn,resi,chain,alt,q))')
    # prepare a list of atom groups
    # each atom group consists of one or more different conformations of the same atom
    alt_atoms_grouped = []
    for atom in stored.alt_atoms:
        match = False
        for group in alt_atoms_grouped:
            if atom[0] == group[0][0] and atom[1] == group[0][1] and atom[2] == group[0][2] and atom[3] == group[0][3]:
                group.append(atom)
                match = True
        if not match:
            alt_atoms_grouped.append([atom])
    # determine which atom conformations to keep, prioritizing the highest occupancy factor and the lowest alt ID, in
    # that order
    atoms_to_remove = []
    for group in alt_atoms_grouped:
        keep = []
        remove = []
        for atom in group:
            if len(keep) == 0:
                keep = atom
            elif atom[5] > keep[5]:
                remove.append(keep)
                keep = atom
            elif atom[5] == keep[5]:
                if len(atom[4]) != 1 or len(keep[4]) != 1:
                    return [False, f"Error: There is at least one atom in equivalence class member {eq_class_mem} with "
                                   f"an alt ID that does not consist of exactly one character."]
                elif atom[4] < keep[4]:
                    remove.append(keep)
                    keep = atom
                elif atom[4] == keep[4]:
                    return [False, f"Error: There are multiple atoms with the same alt ID that should have different "
                                   f"alt IDs in equivalence class member {eq_class_mem}."]
                elif atom[4] > keep[4]:
                    remove.append(atom)
            elif atom[5] < keep[5]:
                remove.append(atom)
        for atom in remove:
            atoms_to_remove.append(atom)
    # remove the atom conformations that are not going to be kept
    for atom in atoms_to_remove:
        # ensure that the info for each atom describes exactly one atom
        if cmd.count_atoms(f'name {atom[0]} and resn {atom[1]} and resi {atom[2]} and chain {atom[3]} and '
                           f'alt {atom[4]}') != 1:
            return [False, f"Error: The info provided for an atom to remove does not account for exactly one atom in "
                           f"equivalence class member {eq_class_mem}."]
        # check for whether there is any indication that removing the atom will disrupt the connectivity of the atoms
        # that are kept
        stored.neighboring_atoms = []
        cmd.iterate(f'neighbor (name {atom[0]} and resn {atom[1]} and resi {atom[2]} and chain {atom[3]} and '
                    f'alt {atom[4]})', 'stored.neighboring_atoms.append((name,resn,resi,chain,alt,q))')
        for neighbor in stored.neighboring_atoms:
            if neighbor[4] != "":
                neighbor_in_remove_list = False
                for compare_atom in atoms_to_remove:
                    if (compare_atom[0] == neighbor[0] and compare_atom[1] == neighbor[1] and
                            compare_atom[2] == neighbor[2] and compare_atom[3] == neighbor[3] and
                            compare_atom[4] == neighbor[4]):
                        neighbor_in_remove_list = True
                if not neighbor_in_remove_list:
                    return [False, f"Error: Removing the atom {atom[3]}.{atom[1]}.{atom[2]}.{atom[0]} with alt ID "
                                   f"{atom[4]} of equivalence class member {eq_class_mem} may disrupt the connectivity "
                                   f"of the atoms that are kept."]
        cmd.remove(f'name {atom[0]} and resn {atom[1]} and resi {atom[2]} and chain {atom[3]} and alt {atom[4]}')
    # Return True as the first element of a list if no errors were encountered.
    return [True]
