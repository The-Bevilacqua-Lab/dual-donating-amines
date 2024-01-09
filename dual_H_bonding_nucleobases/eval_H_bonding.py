"""
This module will eventaully do something.
Mention that PyMOL index is sensitive to atom addition and removal
Mention that histidines should be protonated at both endocyclic nitrogens
"""

import sys
from pymol import cmd
from pymol import stored


# TODO consider replacing the use of index with chain/resi/name
def eval_h_bonding(donor_index, acceptor_index, pdb, h_ang_tol=30, don_ang_tol=60):
    # initially assume that the evaluation will complete successfully
    successful_completion = True
    # initialize an empty list that can be used to document why an evaluation was not completed successfully
    notes = []
    # ensure that the index for the donor atom accounts for exactly one atom
    if cmd.count_atoms(f'index {donor_index}') != 1:
        successful_completion = False
        print(f"Error: The index provided for a potential H-bonding donor atom does not account for exactly one atom "
              f"in PDB ID {pdb}.")
        notes.append(f"Error: The index provided for a potential H-bonding donor atom does not account for exactly one "
                     f"atom in PDB ID {pdb}.")
        return [successful_completion, notes]
    # ensure that the index for the acceptor atom accounts for exactly one atom
    if cmd.count_atoms(f'index {acceptor_index}') != 1:
        successful_completion = False
        print(f"Error: The index provided for a potential H-bonding acceptor atom does not account for exactly one "
              f"atom in PDB ID {pdb}.")
        notes.append(f"Error: The index provided for a potential H-bonding acceptor atom does not account for exactly "
                     f"one atom in PDB ID {pdb}.")
        return [successful_completion, notes]
    # collect information on the donor and acceptor atoms
    stored.donor_atom = ()
    cmd.iterate(f'index {donor_index}', 'stored.donor_atom = (index, name, resn, resi, chain)')
    stored.acceptor_atom = ()
    cmd.iterate(f'index {acceptor_index}', 'stored.acceptor_atom = (index, name, resn, resi, chain)')
    # construct lists containing the names of the canonical protein and nucleic residues
    protein_residues = ['ALA', 'ASP', 'ASN', 'ARG', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
                        'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    nucleic_residues = ['A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT']
    # determine whether the donor atom is at the end of a chain and is capable of donating an H-bond
    terminal_donating_atom = []
    if cmd.count_atoms(f'not elem H and neighbor index {donor_index}') == 1 and (
            stored.donor_atom[1] == "N" or
            stored.donor_atom[1] == "O5'" or
            stored.donor_atom[1] == "O3'"):
        # an amino acid N atom at the end of the chain is capable of donating an H-bond
        if stored.donor_atom[1] == "N" and stored.donor_atom[2] in protein_residues:
            terminal_donating_atom = ["N", "N", True, 0, "CA"]
        else:
            stored.bonded_atom = ''
            cmd.iterate(f'not elem H and neighbor index {donor_index}', 'stored.bonded_atom = name')
            # only an RNA O5' atom at the 5' end of the chain is capable of donating an H-bond
            if (stored.donor_atom[1] == "O5'" and stored.bonded_atom == "C5'" and stored.donor_atom[2] in
                    nucleic_residues):
                terminal_donating_atom = ["O5'", "O", True, 0, "C5'"]
            # only an RNA O3' atom at the 3' end of the chain is capable of donating an H-bond
            if (stored.donor_atom[1] == "O3'" and stored.bonded_atom == "C3'" and stored.donor_atom[2] in
                    nucleic_residues):
                terminal_donating_atom = ["O3'", "O", True, 0, "C3'"]
    # determine whether the donor atom is listed in the residue library
    # if found, collect other information about that atom
    donor_info = []
    donor_res = ''
    donor_atom = ''
    for residue in residue_library:
        if stored.donor_atom[2] == residue['res']:
            donor_res = residue['res']
            for atom in residue['don']:
                if stored.donor_atom[1] == atom[0]:
                    donor_atom = atom[0]
                    donor_info = atom
                    # break the loop since the atom was found
                    break
            if (terminal_donating_atom and (stored.donor_atom[1] == "N" and stored.donor_atom[2] in protein_residues) or
                    ((stored.donor_atom[1] == "O5'" or stored.donor_atom[1] == "O3'") and stored.donor_atom[2] in
                     nucleic_residues)):
                donor_atom = True
                donor_info = terminal_donating_atom
            # break the loop since the residue was found
            break
    if not donor_res:
        successful_completion = False
        notes.append(f"Residue {stored.donor_atom[2]} of the donor atom was not found in the residue library when "
                     f"evaluating a potential H-bond in PDB ID {pdb}.")
    if not donor_atom:
        successful_completion = False
        notes.append(f"The donor atom {stored.donor_atom[1]} of residue {stored.donor_atom[2]} was not found in the "
                     f"residue library when evaluating a potential H-bond in PDB ID {pdb}.")
    # determine whether the acceptor atom is listed in the residue library
    # if found, collect other information about that atom
    acceptor_info = []
    acceptor_res = ''
    acceptor_atom = ''
    for residue in residue_library:
        if stored.acceptor_atom[2] == residue['res']:
            acceptor_res = residue['res']
            for atom in residue['acc']:
                if stored.acceptor_atom[1] == atom[0]:
                    acceptor_atom = atom[0]
                    acceptor_info = atom
                    # break the loop since the atom was found
                    break
            # break the loop since the residue was found
            break
    if not acceptor_res:
        successful_completion = False
        notes.append(f"Residue {stored.acceptor_atom[2]} of the acceptor atom was not found in the residue library "
                     f"when evaluating a potential H-bond in PDB ID {pdb}.")
    if not acceptor_atom:
        successful_completion = False
        notes.append(f"The acceptor atom {stored.acceptor_atom[1]} of residue {stored.acceptor_atom[2]} was not found "
                     f"in the residue library when evaluating a potential H-bond in PDB ID {pdb}.")
    # if any of the residues of the donor or acceptor atoms or the atoms themselves were not found in the residue
    # library, return a non-successful completion with the relevant notes
    if not successful_completion:
        return [successful_completion, notes]
    # if the donor is not rotatable, the vertex of the angle calculation will be the hydrogen atom
    if not donor_info[2]:
        vertex = 'hydrogen'
    # if the donor is rotatable, the vertex of the angle calculation will be the donor atom
    elif donor_info[2]:
        vertex = 'donor'
    # check if the donor and/or acceptor residues have side chains that have conformations that may have been
    # ambiguously assigned
    if donor_res in ["ASN", "GLN", "HIS"] and donor_atom != "N":
        ambiguous_donor = True
    else:
        ambiguous_donor = False
    if acceptor_res in ["ASN", "GLN", "HIS"] and acceptor_atom != "O":
        ambiguous_acceptor = True
    else:
        ambiguous_acceptor = False
    # if there is no ambiguity in the side chain conformations of either donor or acceptor residues, the geometry
    # calculation only needs to be performed once
    if not ambiguous_donor and not ambiguous_acceptor:
        # get the distance and angle values for the donor/acceptor pair
        geometry = calc_geom(donor_index, acceptor_index, donor_info, pdb)
        # if the geometry calculation was not successful, return an explanation of what went wrong
        successful_completion = geometry[0]
        if not successful_completion:
            return geometry
        # proceed if the geometry calculation was successful
        else:
            distance = geometry[1][0]
            angle = geometry[1][1]
            # return whether an H-bond is identified
            if vertex == 'hydrogen':
                # noinspection PyTypeChecker
                if distance <= 2.5 and angle >= (180 - don_ang_tol):
                    return [successful_completion, [True, distance, angle, 'hydrogen']]
                else:
                    return [successful_completion, [False, distance, angle, 'hydrogen']]
            elif vertex == 'donor':
                # noinspection PyTypeChecker
                if distance <= 3.5 and (109.5 - h_ang_tol) <= angle <= (109.5 + h_ang_tol):
                    return [successful_completion, [True, distance, angle, 'donor']]
                else:
                    return [successful_completion, [False, distance, angle, 'donor']]
    # if there is ambiguity in the side chain conformation of either the donor residue or the acceptor residue (but not
    # both), the geometry calculation needs to be performed twice
    elif (ambiguous_donor and not ambiguous_acceptor) or (not ambiguous_donor and ambiguous_acceptor):
        # get the first set of distance and angle values for the donor/acceptor pair
        first_geometry = calc_geom(donor_index, acceptor_index, donor_info, pdb)
        # if the geometry calculation was not successful, return an explanation of what went wrong
        successful_completion = first_geometry[0]
        if not successful_completion:
            return first_geometry
        # proceed if the geometry calculation was successful
        else:
            distance = first_geometry[1][0]
            angle = first_geometry[1][1]
            first_eval = []
            # determine whether an H-bond is identified
            if vertex == 'hydrogen':
                # noinspection PyTypeChecker
                if distance <= 2.5 and angle >= (180 - don_ang_tol):
                    first_eval = [successful_completion, [True, distance, angle, 'hydrogen']]
                else:
                    first_eval = [successful_completion, [False, distance, angle, 'hydrogen']]
            elif vertex == 'donor':
                # noinspection PyTypeChecker
                if distance <= 3.5 and (109.5 - h_ang_tol) <= angle <= (109.5 + h_ang_tol):
                    first_eval = [successful_completion, [True, distance, angle, 'donor']]
                else:
                    first_eval = [successful_completion, [False, distance, angle, 'donor']]
        # rotate the ambiguous side chain atoms of the donor residue
        rotation = rotate_side_chain(stored.donor_atom, pdb)
        # if the rotation was not successful, return an explanation of what went wrong
        successful_completion = rotation[0]
        if not successful_completion:
            return rotation
        # get the second set of distance and angle values for the donor/acceptor pair
        second_geometry = calc_geom(donor_index, acceptor_index, donor_info, pdb)
        # if the geometry calculation was not successful, return an explanation of what went wrong
        successful_completion = second_geometry[0]
        if not successful_completion:
            return second_geometry
        # proceed if the geometry calculation was successful
        else:
            distance = second_geometry[1][0]
            angle = second_geometry[1][1]
            second_eval = []
            # determine whether an H-bond is identified
            if vertex == 'hydrogen':
                # noinspection PyTypeChecker
                if distance <= 2.5 and angle >= (180 - don_ang_tol):
                    second_eval = [successful_completion, [True, distance, angle, 'hydrogen']]
                else:
                    second_eval = [successful_completion, [False, distance, angle, 'hydrogen']]
            elif vertex == 'donor':
                # noinspection PyTypeChecker
                if distance <= 3.5 and (109.5 - h_ang_tol) <= angle <= (109.5 + h_ang_tol):
                    second_eval = [successful_completion, [True, distance, angle, 'donor']]
                else:
                    second_eval = [successful_completion, [False, distance, angle, 'donor']]
        # determine which evaluation to return, giving preference to the first evaluation since that is based on how
        # the structure was originally modeled
        if first_eval[1][0] and second_eval[1][0]:
            return first_eval
        elif first_eval[1][0] and not second_eval[1][0]:
            return first_eval
        elif not first_eval[1][0] and second_eval[1][0]:
            return second_eval
        elif not first_eval[1][0] and not second_eval[1][0]:
            return first_eval
    # if there is ambiguity in the side chain conformations of both donor and acceptor residues, do not calculate any
    # geometries and return False for successful_completion
    # this situation should never be encountered in the current study since, at a minimum, an RNA/DNA residue should
    # always be either the donor or acceptor
    # the code could be expanded to evaluate whether an H-bond is identified in this situation
    elif ambiguous_donor and ambiguous_acceptor:
        successful_completion = False
        print(f"Error: Both the potential H-bonding donor and acceptor atoms belong to the side chains of ASN, GLN, or "
              f"HIS in PDB ID {pdb}. The code is not capable of evaluating whether an H-bond is identified in this "
              f"situation.")
        notes.append(f"Error: Both the potential H-bonding donor and acceptor atoms belong to the side chains of ASN, "
                     f"GLN, or HIS in PDB ID {pdb}. The code is not capable of evaluating whether an H-bond is "
                     f"identified in this situation.")
        return [successful_completion, notes]


def rotate_side_chain(atom, pdb):
    # initially assume that the evaluation will complete successfully
    successful_completion = True
    # initialize an empty list that can be used to document why an evaluation was not completed successfully
    notes = []
    # collect information specific to the residue
    if atom[2] == "ASN":
        origin = cmd.get_coords(f'(byres (name {atom[1]} and resn {atom[2]} and resi {atom[3]} and chain {atom[4]})) '
                                f'and name CG', state=0)[0]
        origin_antecedent = cmd.get_coords(f'(byres (name {atom[1]} and resn {atom[2]} and resi {atom[3]} and chain '
                                           f'{atom[4]})) and name CB', state=0)[0]
        axis = [origin[0] - origin_antecedent[0], origin[1] - origin_antecedent[1], origin[2] - origin_antecedent[2]]
        atoms_to_rotate = (f'(byres (name {atom[1]} and resn {atom[2]} and resi {atom[3]} and chain {atom[4]})) and '
                           f'(name ND2 ((neighbor name ND2) and elem H) name OD1)')
        # ensure that atoms_to_rotate accounts for exactly four atoms
        if cmd.count_atoms(atoms_to_rotate) != 4:
            successful_completion = False
            print(f"Error: Atoms selected from an ASN residue to rotate do not number to exactly four atoms in PDB ID "
                  f"{pdb}.")
            notes.append(f"Error: Atoms selected from an ASN residue to rotate do not number to exactly four atoms in "
                         f"PDB ID {pdb}.")
            return [successful_completion, notes]
    elif atom[2] == "GLN":
        origin = cmd.get_coords(f'(byres (name {atom[1]} and resn {atom[2]} and resi {atom[3]} and chain {atom[4]})) '
                                f'and name CD', state=0)[0]
        origin_antecedent = cmd.get_coords(f'(byres (name {atom[1]} and resn {atom[2]} and resi {atom[3]} and chain '
                                           f'{atom[4]})) and name CG', state=0)[0]
        axis = [origin[0] - origin_antecedent[0], origin[1] - origin_antecedent[1], origin[2] - origin_antecedent[2]]
        atoms_to_rotate = (f'(byres (name {atom[1]} and resn {atom[2]} and resi {atom[3]} and chain {atom[4]})) and '
                           f'(name NE2 ((neighbor name NE2) and elem H) name OE1)')
        # ensure that atoms_to_rotate accounts for exactly four atoms
        if cmd.count_atoms(atoms_to_rotate) != 4:
            successful_completion = False
            print(f"Error: Atoms selected from a GLN residue to rotate do not number to exactly four atoms in PDB ID "
                  f"{pdb}.")
            notes.append(f"Error: Atoms selected from a GLN residue to rotate do not number to exactly four atoms in "
                         f"PDB ID {pdb}.")
            return [successful_completion, notes]
    elif atom[2] == "HIS":
        origin = cmd.get_coords(f'(byres (name {atom[1]} and resn {atom[2]} and resi {atom[3]} and chain {atom[4]})) '
                                f'and name CG', state=0)[0]
        origin_antecedent = cmd.get_coords(f'(byres (name {atom[1]} and resn {atom[2]} and resi {atom[3]} and chain '
                                           f'{atom[4]})) and name CB', state=0)[0]
        axis = [origin[0] - origin_antecedent[0], origin[1] - origin_antecedent[1], origin[2] - origin_antecedent[2]]
        atoms_to_rotate = (f'(byres (name {atom[1]} and resn {atom[2]} and resi {atom[3]} and chain {atom[4]})) and '
                           f'(name ND1 ((neighbor name ND1) and elem H) name CE1 name CD2 name NE2 '
                           f'((neighbor name NE2) and elem H))')
        # ensure that atoms_to_rotate accounts for exactly six atoms
        if cmd.count_atoms(atoms_to_rotate) != 6:
            successful_completion = False
            print(f"Error: Atoms selected from a HIS residue to rotate do not number to exactly six atoms in PDB ID "
                  f"{pdb}.")
            notes.append(f"Error: Atoms selected from a HIS residue to rotate do not number to exactly six atoms in "
                         f"PDB ID {pdb}.")
            return [successful_completion, notes]
    else:
        successful_completion = False
        print(f"Error: The code attempted to rotate a residue in PDB ID {pdb} that is not ASN, GLN, or HIS.")
        notes.append(f"Error: The code attempted to rotate a residue in PDB ID {pdb} that is not ASN, GLN, or HIS.")
        return [successful_completion, notes]
    cmd.rotate(axis, 180, selection=atoms_to_rotate, state=0, camera=0, origin=origin)
    return [successful_completion, notes]


def calc_geom(donor_index, acceptor_index, donor_info, pdb):
    # initially assume that the evaluation will complete successfully
    successful_completion = True
    # initialize an empty list that can be used to document why an evaluation was not completed successfully
    notes = []
    # if the donor is non-rotatable, use the locations of the donor, hydrogen, and acceptor atoms to get the distance
    # and angle measurements, then return the values
    if not donor_info[2]:
        # collect the indices of the hydrogens
        stored.hydrogen = []
        cmd.iterate(f'elem H and neighbor index {donor_index}', 'stored.hydrogen.append(index)')
        # issue an error message if there are no hydrogens
        if len(stored.hydrogen) == 0:
            successful_completion = False
            print(f"Error: A non-rotatable donor atom has zero hydrogens in PDB ID {pdb}.")
            notes.append(f"Error: A non-rotatable donor atom has zero hydrogens in PDB ID {pdb}.")
            return [successful_completion, notes]
        # get the measurements for donors that have one hydrogen
        elif len(stored.hydrogen) == 1:
            distance = cmd.get_distance(f'index {stored.hydrogen[0]}', f'index {acceptor_index}')
            angle = cmd.get_angle(f'index {donor_index}', f'index {stored.hydrogen[0]}', f'index {acceptor_index}')
        # get the measurements for donors that have two hydrogens
        elif len(stored.hydrogen) == 2:
            distance_list = []
            angle_list = []
            for h_index in stored.hydrogen:
                distance_list.append(cmd.get_distance(f'index {h_index}', f'index {acceptor_index}'))
                angle_list.append(cmd.get_angle(f'index {donor_index}', f'index {h_index}', f'index {acceptor_index}'))
            if distance_list[0] < distance_list[1]:
                distance = distance_list[0]
                angle = angle_list[0]
            elif distance_list[0] > distance_list[1]:
                distance = distance_list[1]
                angle = angle_list[1]
            else:
                successful_completion = False
                print(f"Error: The two hydrogens on a non-rotatable donor atom are the same distance from the acceptor "
                      f"atom in PDB ID {pdb}.")
                notes.append(f"Error: The two hydrogens on a non-rotatable donor atom are the same distance from the "
                             f"acceptor atom in PDB ID {pdb}.")
                return [successful_completion, notes]
        # issue an error message if there are more than two hydrogens
        elif len(stored.hydrogen) > 2:
            successful_completion = False
            print(f"Error: A non-rotatable donor atom has more than two hydrogens in PDB ID {pdb}.")
            notes.append(f"Error: A non-rotatable donor atom has more than two hydrogens in PDB ID {pdb}.")
            return [successful_completion, notes]
        # return geometry values
        return [successful_completion, [distance, angle]]
    # if the donor is rotatable, use the locations of the donor antecedent, donor, and acceptor atoms to get the
    # distance and angle measurements, then return the values
    elif donor_info[2]:
        # collect the index of the donor antecedent atom
        stored.antecedent = []
        cmd.iterate(f'name {donor_info[4]} and neighbor index {donor_index}', 'stored.antecedent.append(index)')
        # issue an error message if no donor antecedent atom was identified
        if len(stored.antecedent) == 0:
            successful_completion = False
            print(f"Error: No antecedent atoms were identified for a donor atom in PDB ID {pdb}.")
            notes.append(f"Error: No antecedent atoms were identified for a donor atom in PDB ID {pdb}.")
            return [successful_completion, notes]
        # get the measurements when one donor antecedent atom was identified
        elif len(stored.antecedent) == 1:
            distance = cmd.get_distance(f'index {donor_index}', f'index {acceptor_index}')
            angle = cmd.get_angle(f'index {stored.antecedent[0]}', f'index {donor_index}', f'index {acceptor_index}')
        # issue an error message if more than one donor antecedent atom was identified
        elif len(stored.antecedent) > 1:
            successful_completion = False
            print(f"Error: More than one antecedent atom was identified for a donor atom in PDB ID {pdb}.")
            notes.append(f"Error: More than one antecedent atom was identified for a donor atom in PDB ID {pdb}.")
            return [successful_completion, notes]
        # return geometry values
        return [successful_completion, [distance, angle]]

# region library
# define a list of dictionaries that provides information on the canonical protein and RNA/DNA residues
# the indices of each donor or acceptor atom list specifies the following:
#     0: donor/acceptor atom name
#     1: donor/acceptor element
#     2: is the donor/acceptor rotatable?
#     3: donor/acceptor formal charge
#     4: antecedent atom of donor/acceptor
residue_library = [
    {
        "res": "ALA",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"]]
    },
    {
        "res": "ASP",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"], ["OD1", "O", False, 0, "CG"], ["OD2", "O", False, -1, "CG"]]
    },
    {
        "res": "ASN",
        "don": [["N", "N", False, 0, "CA"], ["ND2", "N", False, 0, "CG"]],
        "acc": [["O", "O", False, 0, "C"], ["OD1", "O", False, 0, "CG"]]
    },
    {
        "res": "ARG",
        "don": [["N", "N", False, 0, "CA"], ["NE", "N", False, 0, "CD"], ["NH1", "N", False, 1, "CZ"],
                ["NH2", "N", False, 0, "CZ"]],
        "acc": [["O", "O", False, 0, "C"]]
    },
    {
        "res": "CYS",
        "don": [["N", "N", False, 0, "CA"], ["SG", "S", True, 0, "CB"]],
        "acc": [["O", "O", False, 0, "C"], ["SG", "S", True, 0, "CB"]]
    },
    {
        "res": "GLU",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"], ["OE1", "O", False, 0, "CD"], ["OE2", "O", False, -1, "CD"]]
    },
    {
        "res": "GLN",
        "don": [["N", "N", False, 0, "CA"], ["NE2", "N", False, 0, "CD"]],
        "acc": [["O", "O", False, 0, "C"], ["OE1", "O", False, 0, "CD"]]
    },
    {
        "res": "GLY",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"]]
    },
    {
        "res": "HIS",
        "don": [["N", "N", False, 0, "CA"], ["ND1", "N", False, 0, "CG"], ["NE2", "N", False, 0, "CD2"]],
        "acc": [["O", "O", False, 0, "C"], ["ND1", "N", False, 0, "CG"], ["NE2", "N", False, 0, "CD2"]]
    },
    {
        "res": "ILE",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"]]
    },
    {
        "res": "LEU",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"]]
    },
    {
        "res": "LYS",
        "don": [["N", "N", False, 0, "CA"], ["NZ", "N", True, 1, "CE"]],
        "acc": [["O", "O", False, 0, "C"]]
    },
    {
        "res": "MET",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"], ["SD", "S", True, 0, "CG"]]
    },
    {
        "res": "PHE",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"]]
    },
    {
        "res": "PRO",
        "don": [],
        "acc": [["O", "O", False, 0, "C"]]
    },
    {
        "res": "SER",
        "don": [["N", "N", False, 0, "CA"], ["OG", "O", True, 0, "CB"]],
        "acc": [["O", "O", False, 0, "C"], ["OG", "O", True, 0, "CB"]]
    },
    {
        "res": "THR",
        "don": [["N", "N", False, 0, "CA"], ["OG1", "O", True, 0, "CB"]],
        "acc": [["O", "O", False, 0, "C"], ["OG1", "O", True, 0, "CB"]]
    },
    {
        "res": "TRP",
        "don": [["N", "N", False, 0, "CA"], ["NE1", "N", False, 0, "CD1"]],
        "acc": [["O", "O", False, 0, "C"]]
    },
    {
        "res": "TYR",
        "don": [["N", "N", False, 0, "CA"], ["OH", "O", False, 0, "CZ"]],
        "acc": [["O", "O", False, 0, "C"], ["OH", "O", False, 0, "CZ"]]
    },
    {
        "res": "VAL",
        "don": [["N", "N", False, 0, "CA"]],
        "acc": [["O", "O", False, 0, "C"]]
    },
    {
        "res": "A",
        "don": [["O2'", "O", True, 0, "C2'"], ["N6", "N", False, 0, "C6"]],
        "acc": [["O2'", "O", True, 0, "C2'"], ["O3'", "O", False, 0, "C3'"], ["O4'", "O", False, 0, "C4'"],
                ["O5'", "O", False, 0, "C5'"], ["OP1", "O", False, 0, "P"], ["OP2", "O", False, -1, "P"],
                ["N1", "N", False, 0, "C2"], ["N3", "N", False, 0, "C2"], ["N7", "N", False, 0, "C5"]]
    },
    {
        "res": "C",
        "don": [["O2'", "O", True, 0, "C2'"], ["N4", "N", False, 0, "C4"]],
        "acc": [["O2'", "O", True, 0, "C2'"], ["O3'", "O", False, 0, "C3'"], ["O4'", "O", False, 0, "C4'"],
                ["O5'", "O", False, 0, "C5'"], ["OP1", "O", False, 0, "P"], ["OP2", "O", False, -1, "P"],
                ["O2", "O", False, 0, "C2"], ["N3", "N", False, 0, "C2"]]
    },
    {
        "res": "G",
        "don": [["O2'", "O", True, 0, "C2'"], ["N1", "N", False, 0, "C2"], ["N2", "N", False, 0, "C2"]],
        "acc": [["O2'", "O", True, 0, "C2'"], ["O3'", "O", False, 0, "C3'"], ["O4'", "O", False, 0, "C4'"],
                ["O5'", "O", False, 0, "C5'"], ["OP1", "O", False, 0, "P"], ["OP2", "O", False, -1, "P"],
                ["N3", "N", False, 0, "C2"], ["O6", "O", False, 0, "C6"], ["N7", "N", False, 0, "C5"]]
    },
    {
        "res": "U",
        "don": [["O2'", "O", True, 0, "C2'"], ["N3", "N", False, 0, "C2"]],
        "acc": [["O2'", "O", True, 0, "C2'"], ["O3'", "O", False, 0, "C3'"], ["O4'", "O", False, 0, "C4'"],
                ["O5'", "O", False, 0, "C5'"], ["OP1", "O", False, 0, "P"], ["OP2", "O", False, -1, "P"],
                ["O2", "O", False, 0, "C2"], ["O4", "O", False, 0, "C4"]]
    },
    {
        "res": "DA",
        "don": [["N6", "N", False, 0, "C6"]],
        "acc": [["O3'", "O", False, 0, "C3'"], ["O4'", "O", False, 0, "C4'"],
                ["O5'", "O", False, 0, "C5'"], ["OP1", "O", False, 0, "P"], ["OP2", "O", False, -1, "P"],
                ["N1", "N", False, 0, "C2"], ["N3", "N", False, 0, "C2"], ["N7", "N", False, 0, "C5"]]
    },
    {
        "res": "DC",
        "don": [["N4", "N", False, 0, "C4"]],
        "acc": [["O3'", "O", False, 0, "C3'"], ["O4'", "O", False, 0, "C4'"],
                ["O5'", "O", False, 0, "C5'"], ["OP1", "O", False, 0, "P"], ["OP2", "O", False, -1, "P"],
                ["O2", "O", False, 0, "C2"], ["N3", "N", False, 0, "C2"]]
    },
    {
        "res": "DG",
        "don": [["N1", "N", False, 0, "C2"], ["N2", "N", False, 0, "C2"]],
        "acc": [["O3'", "O", False, 0, "C3'"], ["O4'", "O", False, 0, "C4'"],
                ["O5'", "O", False, 0, "C5'"], ["OP1", "O", False, 0, "P"], ["OP2", "O", False, -1, "P"],
                ["N3", "N", False, 0, "C2"], ["O6", "O", False, 0, "C6"], ["N7", "N", False, 0, "C5"]]
    },
    {
        "res": "DT",
        "don": [["N3", "N", False, 0, "C2"]],
        "acc": [["O3'", "O", False, 0, "C3'"], ["O4'", "O", False, 0, "C4'"],
                ["O5'", "O", False, 0, "C5'"], ["OP1", "O", False, 0, "P"], ["OP2", "O", False, -1, "P"],
                ["O2", "O", False, 0, "C2"], ["O4", "O", False, 0, "C4"]]
    }
]
# endregion

# print(eval_h_bonding(200457, 22940, '6XU8'))
# rotate_side_chain((61155, 'OD1', 'ASN', '56', 'AG'), '6XU8')
