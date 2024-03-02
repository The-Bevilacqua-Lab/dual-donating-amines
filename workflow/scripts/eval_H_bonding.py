"""
This module contains varoius functions that work together and with PyMOL to calculate H-bonding geometry measurements.
The most notable function is evaluate. This function can be used by other modules or scripts to obtain H-bonding
geometry measurements involving specified donor and acceptor atoms. The first two arguments should describe the donor
and acceptor atoms, respectively. For each atom, a tuple containing the atom index, atom name, residue name, residue
number, and chain ID should be provided in that order. The equivalence class name and a residue library should be
provided as the third and forth arguments. The function assumes 1) that hydrogens have been added to non-rotatable donor
atoms, including both endocyclic nitrogens in histidine residues and 2) that tyrosine hydroxyl hydrogens are
approximately planar with the tyrosine side chain rings. When the potential donor or acceptor atom belongs to the side
chain of an ASN, GLN, or HIS residue, hydrogens must have been added to all the nitrogen atoms on the side chain, even
if the potential H-bond does not involve these nitrogens.
"""

from pymol import cmd
from pymol import stored


def terminal_donor(donor_atom):
    # construct lists containing the names of the canonical protein and nucleic residues
    protein_residues = ['ALA', 'ASP', 'ASN', 'ARG', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
                        'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    nucleic_residues = ['A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT']
    # determine whether the donor atom is at the end of a chain and is capable of donating an H-bond
    terminal_donating_atom = []
    if cmd.count_atoms(f'not elem H and neighbor index {donor_atom[0]}') == 1 and (
            donor_atom[1] == "N" or
            donor_atom[1] == "O5'" or
            donor_atom[1] == "O3'"):
        # an amino acid N atom at the end of the chain is capable of donating an H-bond
        if donor_atom[1] == "N" and donor_atom[2] in protein_residues:
            terminal_donating_atom = ["N", "N", True, 0, "CA"]
        else:
            stored.bonded_atom = ''
            cmd.iterate(f'not elem H and neighbor index {donor_atom[0]}', 'stored.bonded_atom = name')
            # only an RNA O5' atom at the 5' end of the chain is capable of donating an H-bond
            if (donor_atom[1] == "O5'" and stored.bonded_atom == "C5'" and donor_atom[2] in
                    nucleic_residues):
                terminal_donating_atom = ["O5'", "O", True, 0, "C5'"]
            # only an RNA O3' atom at the 3' end of the chain is capable of donating an H-bond
            if (donor_atom[1] == "O3'" and stored.bonded_atom == "C3'" and donor_atom[2] in
                    nucleic_residues):
                # look for a nearby phosphorus just in case PyMOL is not detecting a bond when there really should be
                # one (check out PDB ID 1HMH as an example)
                if cmd.count_atoms(f'elem P within 1.8 of index {donor_atom[0]}') == 0:
                    terminal_donating_atom = ["O3'", "O", True, 0, "C3'"]
    return terminal_donating_atom


def evaluate(donor_atom, acceptor_atom, eq_class, library):
    # initially assume that the evaluation will complete successfully
    successful_completion = True
    # initialize an empty list that can be used to document why an evaluation was not completed successfully
    notes = []
    # construct lists containing the names of the canonical protein and nucleic residues
    protein_residues = ['ALA', 'ASP', 'ASN', 'ARG', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
                        'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    nucleic_residues = ['A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT']
    # determine whether the donor atom is at the end of a chain and is capable of donating an H-bond
    terminal_donating_atom = terminal_donor(donor_atom)
    # determine whether the donor atom is listed in the residue library
    # if found, collect other information about that atom
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
        notes.append(f"Residue {donor_atom[2]} of the donor atom was not found in the residue library when "
                     f"evaluating a potential H-bond in equivalence class {eq_class}.")
    if not donor_info:
        successful_completion = False
        notes.append(f"The donor atom {donor_atom[1]} of residue {donor_atom[2]} was not found in the "
                     f"residue library when evaluating a potential H-bond in equivalence class {eq_class}.")
    # determine whether the acceptor atom is listed in the residue library
    # if found, collect other information about that atom
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
        notes.append(f"Residue {acceptor_atom[2]} of the acceptor atom was not found in the residue library "
                     f"when evaluating a potential H-bond in equivalence class {eq_class}.")
    if not acceptor_info:
        successful_completion = False
        notes.append(f"The acceptor atom {acceptor_atom[1]} of residue {acceptor_atom[2]} was not found "
                     f"in the residue library when evaluating a potential H-bond in equivalence class {eq_class}.")
    # if any of the residues of the donor or acceptor atoms or the atoms themselves were not found in the residue
    # library, return a non-successful completion with the relevant notes
    if not successful_completion:
        return [successful_completion, notes]
    # if the donor is not rotatable, the vertex of the angle calculation will be the hydrogen atom
    if not donor_info[2]:
        vertex = 'hydrogen'
    # if the donor is rotatable, the vertex of the angle calculation will be the donor atom
    else:
        vertex = 'donor'
    # check if the donor and/or acceptor residues have side chains that have conformations that may have been
    # ambiguously assigned
    if donor_res in ["ASN", "GLN", "HIS"] and donor_atom[1] != "N":
        ambiguous_donor = True
    else:
        ambiguous_donor = False
    if acceptor_res in ["ASN", "GLN", "HIS"] and acceptor_atom[1] != "O":
        ambiguous_acceptor = True
    else:
        ambiguous_acceptor = False
    # if there is no ambiguity in the side chain conformations of either donor or acceptor residues and the potential
    # donor is not a tyrosine OH atom, the geometry calculation only needs to be performed once
    if not ambiguous_donor and not ambiguous_acceptor and not (donor_atom[1] == "OH" and donor_atom[2] == "TYR"):
        # get the distance and angle values for the donor/acceptor pair
        geometry = calc_geom(donor_atom, acceptor_atom, donor_info, eq_class)
        # if the geometry calculation was not successful, return an explanation of what went wrong
        successful_completion = geometry[0]
        if not successful_completion:
            return geometry
        # proceed if the geometry calculation was successful
        else:
            # return whether an H-bond is identified
            if vertex == 'hydrogen':
                if type(geometry[1][0]) is float:
                    distance = geometry[1][0]
                    angle = geometry[1][1]
                    h_name = geometry[1][2]
                    return [successful_completion, [[distance, angle, 'hydrogen', h_name, 'none']]]
                elif type(geometry[1][0]) is list:
                    first_distance = geometry[1][0][0]
                    first_angle = geometry[1][0][1]
                    first_h_name = geometry[1][0][2]
                    first_set = [first_distance, first_angle, 'hydrogen', first_h_name, 'none']
                    second_distance = geometry[1][1][0]
                    second_angle = geometry[1][1][1]
                    second_h_name = geometry[1][1][2]
                    second_set = [second_distance, second_angle, 'hydrogen', second_h_name, 'none']
                    return [successful_completion, [first_set, second_set]]
            elif vertex == 'donor':
                distance = geometry[1][0]
                angle = geometry[1][1]
                h_name = geometry[1][2]
                return [successful_completion, [[distance, angle, 'donor', h_name, 'none']]]
    # if there is ambiguity in the side chain conformation of either the donor residue or the acceptor residue (but not
    # both), the geometry calculation needs to be performed twice
    elif (((ambiguous_donor and not ambiguous_acceptor) or (not ambiguous_donor and ambiguous_acceptor)) and not
            (donor_atom[1] == "OH" and donor_atom[2] == "TYR")):
        # get the first set of distance and angle values for the donor/acceptor pair
        first_geometry = calc_geom(donor_atom, acceptor_atom, donor_info, eq_class)
        # if the geometry calculation was not successful, return an explanation of what went wrong
        successful_completion = first_geometry[0]
        if not successful_completion:
            return first_geometry
        # proceed if the geometry calculation was successful
        else:
            # determine whether an H-bond is identified
            first_eval = []
            if vertex == 'hydrogen':
                if type(first_geometry[1][0]) is float:
                    distance = first_geometry[1][0]
                    angle = first_geometry[1][1]
                    h_name = first_geometry[1][2]
                    first_eval = [[distance, angle, 'hydrogen', h_name, 'none']]
                elif type(first_geometry[1][0]) is list:
                    first_distance = first_geometry[1][0][0]
                    first_angle = first_geometry[1][0][1]
                    first_h_name = first_geometry[1][0][2]
                    first_set = [first_distance, first_angle, 'hydrogen', first_h_name, 'none']
                    second_distance = first_geometry[1][1][0]
                    second_angle = first_geometry[1][1][1]
                    second_h_name = first_geometry[1][1][2]
                    second_set = [second_distance, second_angle, 'hydrogen', second_h_name, 'none']
                    first_eval = [first_set, second_set]
            elif vertex == 'donor':
                distance = first_geometry[1][0]
                angle = first_geometry[1][1]
                h_name = first_geometry[1][2]
                first_eval = [[distance, angle, 'donor', h_name, 'none']]
        # rotate the ambiguous side chain atoms of the donor or acceptor residue
        if ambiguous_donor:
            rotated_res = f'{donor_atom[2]}{donor_atom[3]} {donor_atom[4]}'
            rotation = rotate_side_chain(donor_atom, eq_class)
        else:
            rotated_res = f'{acceptor_atom[2]}{acceptor_atom[3]} {acceptor_atom[4]}'
            rotation = rotate_side_chain(acceptor_atom, eq_class)
        # if the rotation was not successful, return an explanation of what went wrong
        successful_completion = rotation[0]
        if not successful_completion:
            return rotation
        # get the second set of distance and angle values for the donor/acceptor pair
        second_geometry = calc_geom(donor_atom, acceptor_atom, donor_info, eq_class)
        # if the geometry calculation was not successful, return an explanation of what went wrong
        successful_completion = second_geometry[0]
        if not successful_completion:
            return second_geometry
        # proceed if the geometry calculation was successful
        else:
            # determine whether an H-bond is identified
            second_eval = []
            if vertex == 'hydrogen':
                if type(second_geometry[1][0]) is float:
                    distance = second_geometry[1][0]
                    angle = second_geometry[1][1]
                    h_name = second_geometry[1][2]
                    second_eval = [[distance, angle, 'hydrogen', h_name, rotated_res]]
                elif type(second_geometry[1][0]) is list:
                    first_distance = second_geometry[1][0][0]
                    first_angle = second_geometry[1][0][1]
                    first_h_name = second_geometry[1][0][2]
                    first_set = [first_distance, first_angle, 'hydrogen', first_h_name, rotated_res]
                    second_distance = second_geometry[1][1][0]
                    second_angle = second_geometry[1][1][1]
                    second_h_name = second_geometry[1][1][2]
                    second_set = [second_distance, second_angle, 'hydrogen', second_h_name, rotated_res]
                    second_eval = [first_set, second_set]
            elif vertex == 'donor':
                distance = second_geometry[1][0]
                angle = second_geometry[1][1]
                h_name = second_geometry[1][2]
                second_eval = [[distance, angle, 'donor', h_name, rotated_res]]
        # rotate the ambiguous side chain atoms of the donor or acceptor residue to get them back to their original
        # conformation
        if ambiguous_donor:
            rotation = rotate_side_chain(donor_atom, eq_class)
        elif ambiguous_acceptor:
            rotation = rotate_side_chain(acceptor_atom, eq_class)
        # if the rotation was not successful, return an explanation of what went wrong
        successful_completion = rotation[0]
        if not successful_completion:
            return rotation
        # return all evaluations
        all_eval = first_eval
        all_eval.extend(second_eval)
        return [successful_completion, all_eval]
    # if there is ambiguity in the side chain conformations of both donor and acceptor residues, do not calculate any
    # geometries and return False for successful_completion
    # this situation should never be encountered in the current study since, at a minimum, an RNA/DNA residue should
    # always be either the donor or acceptor
    # the code could be expanded to evaluate whether an H-bond is identified in this situation
    elif ambiguous_donor and ambiguous_acceptor and not (donor_atom[1] == "OH" and donor_atom[2] == "TYR"):
        successful_completion = False
        notes.append(f"Error: Both the potential H-bonding donor and acceptor atoms belong to the side chains of ASN, "
                     f"GLN, or HIS in equivalence class {eq_class}. The code is not capable of evaluating whether an "
                     f"H-bond is identified in this situation.")
        return [successful_completion, notes]
    elif donor_atom[1] == "OH" and donor_atom[2] == "TYR":
        # get the distance and angle values for the donor/acceptor pair and with the tyrosine OH in the original
        # conformation
        geometry = calc_geom(donor_atom, acceptor_atom, donor_info, eq_class)
        # if the geometry calculation was not successful, return an explanation of what went wrong
        successful_completion = geometry[0]
        if not successful_completion:
            return geometry
        # determine whether an H-bond is identified
        distance = geometry[1][0]
        angle = geometry[1][1]
        h_name = geometry[1][2]
        original_conformation = [distance, angle, 'hydrogen', h_name, 'none']
        # rotate the tyrosine hydroxyl
        rotated_res = f'{donor_atom[2]}{donor_atom[3]} {donor_atom[4]}'
        rotation = rotate_side_chain(donor_atom, eq_class)
        # if the rotation was not successful, return an explanation of what went wrong
        successful_completion = rotation[0]
        if not successful_completion:
            return rotation
        # get the distance and angle values for the donor/acceptor pair and with the tyrosine OH in the rotated
        # conformation
        geometry = calc_geom(donor_atom, acceptor_atom, donor_info, eq_class)
        # if the geometry calculation was not successful, return an explanation of what went wrong
        successful_completion = geometry[0]
        if not successful_completion:
            return geometry
        # determine whether an H-bond is identified
        distance = geometry[1][0]
        angle = geometry[1][1]
        h_name = geometry[1][2]
        rotated_conformation = [distance, angle, 'hydrogen', h_name, rotated_res]
        # rotate the tyrosine hydroxyl to get it back to its original conformation
        rotation = rotate_side_chain(donor_atom, eq_class)
        # if the rotation was not successful, return an explanation of what went wrong
        successful_completion = rotation[0]
        if not successful_completion:
            return rotation
        # return all evaluations
        return [successful_completion, [original_conformation, rotated_conformation]]
    # I don't know why the code would ever get to this point, but this will serve as a catch-all for anything that I may
    # have missed
    else:
        successful_completion = False
        notes.append(f"Error: The code was unable to classify the donor and/or acceptor as non-ambiguous, ambiguous, "
                     f"or Tyr(OH).")
        return [successful_completion, notes]


def rotate_side_chain(atom, eq_class):
    # initially assume that the evaluation will complete successfully
    successful_completion = True
    # initialize an empty list that can be used to document why an evaluation was not completed successfully
    notes = []
    # collect information specific to the residue
    if atom[2] == "ASN":
        origin = cmd.get_coords(f'(byres index {atom[0]}) and name CG', state=0)[0]
        origin_antecedent = cmd.get_coords(f'(byres index {atom[0]}) and name CB', state=0)[0]
        axis = [origin[0] - origin_antecedent[0], origin[1] - origin_antecedent[1], origin[2] - origin_antecedent[2]]
        atoms_to_rotate = f'(byres index {atom[0]}) and (name ND2 ((neighbor name ND2) and elem H) name OD1)'
        # ensure that atoms_to_rotate accounts for exactly four atoms
        if cmd.count_atoms(atoms_to_rotate) != 4:
            successful_completion = False
            notes.append(f"Error: Atoms selected from an ASN residue to rotate do not number to exactly four atoms in "
                         f"equivalence class {eq_class}.")
            return [successful_completion, notes]
    elif atom[2] == "GLN":
        origin = cmd.get_coords(f'(byres index {atom[0]}) and name CD', state=0)[0]
        origin_antecedent = cmd.get_coords(f'(byres index {atom[0]}) and name CG', state=0)[0]
        axis = [origin[0] - origin_antecedent[0], origin[1] - origin_antecedent[1], origin[2] - origin_antecedent[2]]
        atoms_to_rotate = f'(byres index {atom[0]}) and (name NE2 ((neighbor name NE2) and elem H) name OE1)'
        # ensure that atoms_to_rotate accounts for exactly four atoms
        if cmd.count_atoms(atoms_to_rotate) != 4:
            successful_completion = False
            notes.append(f"Error: Atoms selected from a GLN residue to rotate do not number to exactly four atoms in "
                         f"equivalence class {eq_class}.")
            return [successful_completion, notes]
    elif atom[2] == "HIS":
        origin = cmd.get_coords(f'(byres index {atom[0]}) and name CG', state=0)[0]
        origin_antecedent = cmd.get_coords(f'(byres index {atom[0]}) and name CB', state=0)[0]
        axis = [origin[0] - origin_antecedent[0], origin[1] - origin_antecedent[1], origin[2] - origin_antecedent[2]]
        atoms_to_rotate = (f'(byres index {atom[0]}) and (name ND1 ((neighbor name ND1) and elem H) name CE1 name CD2 '
                           f'name NE2 ((neighbor name NE2) and elem H))')
        # ensure that atoms_to_rotate accounts for exactly six atoms
        if cmd.count_atoms(atoms_to_rotate) != 6:
            successful_completion = False
            notes.append(f"Error: Atoms selected from a HIS residue to rotate do not number to exactly six atoms in "
                         f"equivalence class {eq_class}.")
            return [successful_completion, notes]
    elif atom[2] == "TYR":
        origin = cmd.get_coords(f'index {atom[0]}', state=0)[0]
        origin_antecedent = cmd.get_coords(f'(neighbor index {atom[0]}) and name CZ', state=0)[0]
        axis = [origin[0] - origin_antecedent[0], origin[1] - origin_antecedent[1], origin[2] - origin_antecedent[2]]
        atoms_to_rotate = f'index {atom[0]} ((neighbor index {atom[0]}) and elem H)'
        # ensure that atoms_to_rotate accounts for exactly two atoms
        if cmd.count_atoms(atoms_to_rotate) != 2:
            successful_completion = False
            notes.append(f"Error: Atoms selected from a TYR residue to rotate do not number to exactly two atoms in "
                         f"equivalence class {eq_class}.")
            return [successful_completion, notes]
    else:
        successful_completion = False
        notes.append(f"Error: The code attempted to rotate a residue in equivalence class {eq_class} that is not ASN, "
                     f"GLN, or HIS.")
        return [successful_completion, notes]
    cmd.rotate(axis, 180, selection=atoms_to_rotate, state=0, camera=0, origin=origin)
    return [successful_completion, notes]


def calc_geom(donor_atom, acceptor_atom, donor_info, eq_class):
    # initially assume that the evaluation will complete successfully
    successful_completion = True
    # initialize an empty list that can be used to document why an evaluation was not completed successfully
    notes = []
    # if the donor is non-rotatable, use the locations of the donor, hydrogen, and acceptor atoms to get the distance
    # and angle measurements, then return the values
    if not donor_info[2]:
        # collect the indices of the hydrogens
        stored.hydrogen = []
        cmd.iterate(f'elem H and neighbor index {donor_atom[0]}', 'stored.hydrogen.append((index, name, resn, resi, '
                                                                  'chain))')
        # issue an error message if there are no hydrogens
        if len(stored.hydrogen) == 0:
            successful_completion = False
            notes.append(f"Error: A non-rotatable donor atom has zero hydrogens in equivalence class {eq_class}.")
            return [successful_completion, notes]
        # get the measurements for donors that have one hydrogen
        elif len(stored.hydrogen) == 1:
            distance = cmd.get_distance(f'index {stored.hydrogen[0][0]}', f'index {acceptor_atom[0]}')
            angle = cmd.get_angle(f'index {donor_atom[0]}', f'index {stored.hydrogen[0][0]}',
                                  f'index {acceptor_atom[0]}')
            # return geometry values
            return [successful_completion, [distance, angle, stored.hydrogen[0][1]]]
        # get the measurements for donors that have two hydrogens
        elif len(stored.hydrogen) == 2:
            first_distance = cmd.get_distance(f'index {stored.hydrogen[0][0]}', f'index {acceptor_atom[0]}')
            first_angle = cmd.get_angle(f'index {donor_atom[0]}', f'index {stored.hydrogen[0][0]}',
                                        f'index {acceptor_atom[0]}')
            second_distance = cmd.get_distance(f'index {stored.hydrogen[1][0]}', f'index {acceptor_atom[0]}')
            second_angle = cmd.get_angle(f'index {donor_atom[0]}', f'index {stored.hydrogen[1][0]}',
                                         f'index {acceptor_atom[0]}')
            # return both sets of geometry values
            return [successful_completion, [[first_distance, first_angle, stored.hydrogen[0][1]],
                                            [second_distance, second_angle, stored.hydrogen[1][1]]]]
        # issue an error message if there are more than two hydrogens
        elif len(stored.hydrogen) > 2:
            successful_completion = False
            notes.append(f"Error: A non-rotatable donor atom has more than two hydrogens in equivalence class "
                         f"{eq_class}.")
            return [successful_completion, notes]
    # if the donor is rotatable, use the locations of the donor antecedent, donor, and acceptor atoms to get the
    # distance and angle measurements, then return the values
    elif donor_info[2]:
        # collect information about the donor antecedent atom
        stored.antecedent = []
        cmd.iterate(f'name {donor_info[4]} and neighbor index {donor_atom[0]}',
                    'stored.antecedent.append((index, name, resn, resi, chain))')
        # issue an error message if no donor antecedent atom was identified
        if len(stored.antecedent) == 0:
            successful_completion = False
            notes.append(f"Error: No antecedent atoms were identified for a donor atom in equivalence class "
                         f"{eq_class}.")
            return [successful_completion, notes]
        # get the measurements when one donor antecedent atom was identified
        elif len(stored.antecedent) == 1:
            distance = cmd.get_distance(f'index {donor_atom[0]}', f'index {acceptor_atom[0]}')
            angle = cmd.get_angle(f'index {stored.antecedent[0][0]}', f'index {donor_atom[0]}',
                                  f'index {acceptor_atom[0]}')
        # issue an error message if more than one donor antecedent atom was identified
        else:
            successful_completion = False
            notes.append(f"Error: More than one antecedent atom was identified for a donor atom in equivalence class "
                         f"{eq_class}.")
            return [successful_completion, notes]
        # return geometry values
        return [successful_completion, [distance, angle, "none"]]
