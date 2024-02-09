"""
This module reads information on the equivalence classes in the representative set file. For each equivalence class, it
works with PyMOL to save a mmCIF file of the original structure to the original_mmCIF_files folder and another mmCIF
file with modifications made (alternative conformations removed and hydrogens added) to the modified_mmCIF_files folder.
"""

import sys
import os
import csv
from pymol import cmd
from pymol import stored
from datetime import datetime
import parse_nrlist
import residue_library
import eval_H_bonding
import remove_alt_conf

start = datetime.now()

directory = os.getcwd()
original_mmCIF_directory = directory + "/original_mmCIF_files"
modified_mmCIF_directory = directory + "/modified_mmCIF_files"

# identify the representative set file
nrlist_file = parse_nrlist.identify()

# iterate through the lines of the representative set file and collect the equivalence class names and PDB IDs, model
# info, and chain info for the representative structures
nrlist_info = parse_nrlist.get_info(nrlist_file)

# # check whether an original_mmCIF_files folder exists
# # if it exists and contains files, exit with an error message
# # if it does not exist, create the folder
# if os.path.isdir(original_mmCIF_directory):
#     if len(os.listdir(original_mmCIF_directory)) != 0:
#         print("Error: There are already files in the original_mmCIF_files folder. This folder should be empty prior to "
#               "running the prepare_PDB module.")
#         sys.exit(1)
# else:
#     os.mkdir(original_mmCIF_directory)
#
# # change fetch_path to the original_mmCIF_files folder in the current working directory so that PyMOL drops the mmCIF
# # files into this folder when running fetch
# cmd.set('fetch_path', cmd.exp_path('./original_mmCIF_files'))
#
# # check whether a modified_mmCIF_files folder exists
# # if it exists and contains files, exit with an error message
# # if it does not exist, create the folder
# if os.path.isdir(modified_mmCIF_directory):
#     if len(os.listdir(modified_mmCIF_directory)) != 0:
#         print("Error: There are already files in the modified_mmCIF_files folder. This folder should be empty prior to "
#               "running the prepare_PDB module.")
#         sys.exit(1)
# else:
#     os.mkdir(modified_mmCIF_directory)

# construct two lists of tuples containing atom and residue names that describe either donor or acceptor atoms
# additionally, construct a list of tuples containing atom and residue names that describe rotatable donor atoms
donor_atoms = []
rotatable_donor_atoms = []
acceptor_atoms = []
for residue in residue_library.residue_library:
    for donor in residue['don']:
        donor_atoms.append((residue['res'], donor[0]))
        if donor[2]:
            rotatable_donor_atoms.append((residue['res'], donor[0]))
    for acceptor in residue['acc']:
        acceptor_atoms.append((residue['res'], acceptor[0]))

# construct three strings that can be used with PyMOL to select all possible donor, all possible rotatable donor, or
# all possible acceptor atoms
donor_string = ''
last_atom = donor_atoms[-1]
for atom in donor_atoms:
    donor_string += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        donor_string += ' '
rotatable_donor_string = ''
last_atom = rotatable_donor_atoms[-1]
for atom in rotatable_donor_atoms:
    rotatable_donor_string += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        rotatable_donor_string += ' '
acceptor_string = ''
last_atom = acceptor_atoms[-1]
for atom in acceptor_atoms:
    acceptor_string += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        acceptor_string += ' '

# construct two lists of tuples containing atom and residue names that describe either donor or acceptor atoms that
# belong to side chains that have conformations that may have been ambiguously assigned
ambiguous_donor_atoms = []
ambiguous_acceptor_atoms = []
for residue in residue_library.residue_library:
    if residue['res'] in ["ASN", "GLN", "HIS"]:
        for donor in residue['don']:
            if donor[0] != "N":
                ambiguous_donor_atoms.append((residue['res'], donor[0]))
        for acceptor in residue['acc']:
            if acceptor[0] != "O":
                ambiguous_acceptor_atoms.append((residue['res'], acceptor[0]))

# construct two strings that can be used with PyMOL to select all possible donor or all possible acceptor atoms that
# belong to side chains that have conformations that may have been ambiguously assigned
ambiguous_donor_string = ''
last_atom = ambiguous_donor_atoms[-1]
for atom in ambiguous_donor_atoms:
    ambiguous_donor_string += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        ambiguous_donor_string += ' '
ambiguous_acceptor_string = ''
last_atom = ambiguous_acceptor_atoms[-1]
for atom in ambiguous_acceptor_atoms:
    ambiguous_acceptor_string += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        ambiguous_acceptor_string += ' '

# construct three tuples of tuples containing atom and residue names that describe either donor, protonated donor
# (formal charge of +1), or acceptor atoms of particular interest
donors_of_interest = (('A', 'N6'), ('C', 'N4'), ('G', 'N2'))
prot_donors_of_interest = (('A', 'N1'), ('A', 'N3'), ('C', 'N3'), ('G', 'N3'))
acceptors_of_interest = (('A', 'N1'), ('A', 'N3'), ('C', 'O2'), ('C', 'N3'), ('G', 'N3'), ('G', 'O6'))

# construct three strings that can be used with PyMOL to select all possible donor, all possible protonated donor
# (formal charge of +1), or all possible acceptor atoms of particular interest
donors_of_interest_str = ''
last_atom = donors_of_interest[-1]
for atom in donors_of_interest:
    donors_of_interest_str += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        donors_of_interest_str += ' '
prot_donors_of_interest_str = ''
last_atom = prot_donors_of_interest[-1]
for atom in prot_donors_of_interest:
    prot_donors_of_interest_str += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        prot_donors_of_interest_str += ' '
acceptors_of_interest_str = ''
last_atom = acceptors_of_interest[-1]
for atom in acceptors_of_interest:
    acceptors_of_interest_str += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        acceptors_of_interest_str += ' '

# define the donor and acceptor atoms of the RNA nucleobases
adenine_nuc_donors = ["N6"]
adenine_nuc_acceptors = ["N1", "N3", "N7"]
cytosine_nuc_donors = ["N4"]
cytosine_nuc_acceptors = ["O2", "N3"]
guanine_nuc_donors = ["N1", "N2"]
guanine_nuc_acceptors = ["N3", "O6", "N7"]
uracil_nuc_donors = ["N3"]
uracil_nuc_acceptors = ["O2", "O4"]

# using the library provided by the residue_library module, create a new library with protonated A, C, and G residues
prot_library = residue_library.residue_library
for res in prot_library:
    if res['res'] == 'A':
        res['don'].append(["N1", "N", False, 1, "C2"])
        res['don'].append(["N3", "N", False, 1, "C2"])
    if res['res'] == 'C':
        res['don'].append(["N3", "N", False, 1, "C2"])
    if res['res'] == 'G':
        res['don'].append(["N3", "N", False, 1, "C2"])

# define the distance used to search for nearby donor or acceptor atoms
search_dist = 3.6
# define the distance used to search for nearby donor or acceptor atoms that belong to side chains that have
# conformations that may have been ambiguously assigned
search_dist_amb = search_dist + 3.0

# TODO build in some sort of check to insure that index values do not change when they are expected to remain constant
# work with each equivalence class specified in the representative set file
for eq_class in nrlist_info:
    rep_rna_list = []
    for i in range(len(eq_class[2])):
        # the create function below assumes that the different representative RNA chains belong to the same model
        # if this is not the case, print an error message and exit
        if eq_class[2][i] != eq_class[2][0]:
            print(f"Error: There are representative RNA chains that belong to different models for equivalence class "
                  f"{eq_class[0]}.")
            sys.exit(1)
        rep_rna_list.append(f"chain {eq_class[3][i]}")
    rep_rna = " ".join(rep_rna_list)
    # # delete all objects in the current PyMOL session
    # cmd.delete('all')
    # # retrieve the structure that contains the representative RNA chains
    # cmd.fetch(eq_class[1][0])
    # if the structure contains more than one model (or "state" in PyMOL), create a new PyMOL object of just that model
    if cmd.count_states(eq_class[1][0]) > 1:
        cmd.create(f'{eq_class[1][0]}_state_{eq_class[2][0]}', selection=eq_class[1][0], source_state=eq_class[2][0],
                   target_state=1)
        cmd.delete(eq_class[1][0])
    # remove atoms representing alternative conformations
    remove_status = remove_alt_conf.remove(eq_class[1][0])
    successful_completion = remove_status[0]
    if not successful_completion:
        for note in remove_status[1]:
            print(note)
        sys.exit(1)

    # The h_add command a few lines down will instruct PyMOL to add hydrogens to tyrosine OH atoms (and other atoms
    # too). The tyrosine OH atoms need to be trigonal planar and the hydrogens added to the tyrosine OH atoms need to be
    # planar with the tyrosine rings. For this to happen, there needs to be a double bond between the tyrosine OH and CZ
    # atoms and the tyrosine OH atom needs to have a formal charge of +1.

    # create a double bond between the CZ and OH atoms in tyrosine residues that are near the representative RNA
    cmd.iterate(f'(resn TYR and name OH) within {search_dist_amb} of ({rep_rna})',
                'cmd.unbond(f"index {index}", f"(byres index {index}) and name CZ"); '
                'cmd.bond(f"index {index}", f"(byres index {index}) and name CZ", "2")')
    # change the formal charge on TYR(OH), HIS(ND1), A(N1), A(N3), C(N3), and G(N3) to +1
    # with a formal charge of +1, PyMOL will add hydrogens to these atoms when the following h_add command is used
    cmd.alter(f'(resn TYR and name OH) (resn HIS and name ND1) {prot_donors_of_interest_str}', 'formal_charge=1')
    # add hydrogens to non-rotatable donors that are a part of or near the representative RNA
    cmd.h_add(f'(({donor_string} {prot_donors_of_interest_str}) and not ({rotatable_donor_string})) within '
              f'{search_dist_amb + 3.0} of ({rep_rna})')
    # store a list of donors of interest from the representative RNA
    stored.donor_list = []
    cmd.iterate(f'({rep_rna}) and ({donors_of_interest_str})',
                'stored.donor_list.append((index, name, resn, resi, chain))')
    # store a list of atoms near the donors of interest
    atoms_near_donors = []
    # stored.donor_list = stored.donor_list[549:580]
    for donor in stored.donor_list:
        stored.nearby_atoms = []
        cmd.iterate(f'index {donor[0]} around {search_dist}',
                    'stored.nearby_atoms.append((index, name, resn, resi, chain))')
        atoms_near_donors.append(stored.nearby_atoms)
    # extract the atoms that can act as H-bond acceptors
    acceptors_near_donors = []
    for atom_group in enumerate(atoms_near_donors):
        list_of_acceptors = []
        for atom in atom_group[1]:
            # the side chain acceptor atoms in ASN, GLN, and HIS residues will be collected separately
            if not (atom[1] != "O" and atom[2] in ["ASN", "GLN", "HIS"]):
                # the acceptors should not belong to the nucleobase of the donor
                if not ((stored.donor_list[atom_group[0]][2] == atom[2] == "A" and
                         stored.donor_list[atom_group[0]][3] == atom[3] and
                         stored.donor_list[atom_group[0]][4] == atom[4] and
                         atom[1] in adenine_nuc_acceptors) or
                        (stored.donor_list[atom_group[0]][2] == atom[2] == "C" and
                         stored.donor_list[atom_group[0]][3] == atom[3] and
                         stored.donor_list[atom_group[0]][4] == atom[4] and
                         atom[1] in cytosine_nuc_acceptors) or
                        (stored.donor_list[atom_group[0]][2] == atom[2] == "G" and
                         stored.donor_list[atom_group[0]][3] == atom[3] and
                         stored.donor_list[atom_group[0]][4] == atom[4] and
                         atom[1] in guanine_nuc_acceptors) or
                        (stored.donor_list[atom_group[0]][2] == atom[2] == "U" and
                         stored.donor_list[atom_group[0]][3] == atom[3] and
                         stored.donor_list[atom_group[0]][4] == atom[4] and
                         atom[1] in uracil_nuc_acceptors)):
                    for acceptor in acceptor_atoms:
                        if atom[1] == acceptor[1] and atom[2] == acceptor[0]:
                            list_of_acceptors.append(atom)
        acceptors_near_donors.append(list_of_acceptors)
    # acquire the H-bonding geometry measurements for all acceptors near each donor
    donor_h_bonds = []
    for donor in enumerate(stored.donor_list):
        h_bonds = []
        for acceptor in acceptors_near_donors[donor[0]]:
            h_bonds.append(eval_H_bonding.evaluate(donor[1], acceptor, eq_class[1][0], prot_library))
            # if the H-bond evaluation is not successful, print the error message and exit
            successful_completion = h_bonds[-1][0]
            if not successful_completion:
                for note in h_bonds[-1][1]:
                    print(note)
                sys.exit(1)
        donor_h_bonds.append(h_bonds)
    # using a greater search distance for the side chains of ASN, GLN, and HIS residues, store a list of atoms near the
    # donors of interest
    atoms_near_donors_amb = []
    for donor in stored.donor_list:
        stored.nearby_atoms = []
        cmd.iterate(f'index {donor[0]} around {search_dist_amb}',
                    'stored.nearby_atoms.append((index, name, resn, resi, chain))')
        atoms_near_donors_amb.append(stored.nearby_atoms)
    # extract the atoms that can act as H-bond acceptors that belong to the side chains of ASN, GLN, and HIS residues
    acceptors_near_donors_amb = []
    for atom_group in enumerate(atoms_near_donors_amb):
        list_of_acceptors = []
        for atom in atom_group[1]:
            # only extract side chain acceptor atoms in ASN, GLN, and HIS residues
            if atom[1] != "O" and atom[2] in ["ASN", "GLN", "HIS"]:
                for acceptor in acceptor_atoms:
                    if atom[1] == acceptor[1] and atom[2] == acceptor[0]:
                        list_of_acceptors.append(atom)
        acceptors_near_donors_amb.append(list_of_acceptors)
    # acquire the H-bonding geometry measurements for side chain acceptor atoms in ASN, GLN, and HIS residues near each
    # donor
    donor_h_bonds_amb = []
    for donor in enumerate(stored.donor_list):
        h_bonds = []
        for acceptor in acceptors_near_donors_amb[donor[0]]:
            h_bonds.append(eval_H_bonding.evaluate(donor[1], acceptor, eq_class[1][0], prot_library))
            # if the H-bond evaluation is not successful, print the error message and exit
            successful_completion = h_bonds[-1][0]
            if not successful_completion:
                for note in h_bonds[-1][1]:
                    print(note)
                sys.exit(1)
        donor_h_bonds_amb.append(h_bonds)
    # store a list of acceptors of interest from the representative RNA
    stored.acceptor_list = []
    cmd.iterate(f'({rep_rna}) and ({acceptors_of_interest_str})',
                'stored.acceptor_list.append((index, name, resn, resi, chain))')
    # store a list of atoms near the acceptors of interest
    # additionally, use the list of acceptors of interest to construct a list of protonated donors of interest (formal
    # charge of +1) and a list of atoms near the protonated donors of interest
    atoms_near_acceptors = []
    prot_donor_list = []
    atoms_near_prot_donors = []
    # stored.acceptor_list = stored.acceptor_list[:100]
    for acceptor in stored.acceptor_list:
        stored.nearby_atoms = []
        cmd.iterate(f'index {acceptor[0]} around {search_dist}',
                    'stored.nearby_atoms.append((index, name, resn, resi, chain))')
        atoms_near_acceptors.append(stored.nearby_atoms)
        if (acceptor[2], acceptor[1]) in prot_donors_of_interest:
            prot_donor_list.append(acceptor)
            atoms_near_prot_donors.append(stored.nearby_atoms)
    # extract the atoms near the acceptors of interest that can act as H-bond donors
    donors_near_acceptors = []
    for atom_group in enumerate(atoms_near_acceptors):
        list_of_donors = []
        for atom in atom_group[1]:
            # the side chain donor atoms in ASN, GLN, and HIS residues will be collected separately
            if not (atom[1] != "N" and atom[2] in ["ASN", "GLN", "HIS"]):
                # the donors should not belong to the nucleobase of the acceptor
                if not ((stored.acceptor_list[atom_group[0]][2] == atom[2] == "A" and
                         stored.acceptor_list[atom_group[0]][3] == atom[3] and
                         stored.acceptor_list[atom_group[0]][4] == atom[4] and
                         atom[1] in adenine_nuc_donors) or
                        (stored.acceptor_list[atom_group[0]][2] == atom[2] == "C" and
                         stored.acceptor_list[atom_group[0]][3] == atom[3] and
                         stored.acceptor_list[atom_group[0]][4] == atom[4] and
                         atom[1] in cytosine_nuc_donors) or
                        (stored.acceptor_list[atom_group[0]][2] == atom[2] == "G" and
                         stored.acceptor_list[atom_group[0]][3] == atom[3] and
                         stored.acceptor_list[atom_group[0]][4] == atom[4] and
                         atom[1] in guanine_nuc_donors) or
                        (stored.acceptor_list[atom_group[0]][2] == atom[2] == "U" and
                         stored.acceptor_list[atom_group[0]][3] == atom[3] and
                         stored.acceptor_list[atom_group[0]][4] == atom[4] and
                         atom[1] in uracil_nuc_donors)):
                    for donor in donor_atoms:
                        if atom[1] == donor[1] and atom[2] == donor[0]:
                            list_of_donors.append(atom)
        donors_near_acceptors.append(list_of_donors)
    # acquire the H-bonding geometry measurements for all donors near each acceptor
    acceptor_h_bonds = []
    for acceptor in enumerate(stored.acceptor_list):
        h_bonds = []
        for donor in donors_near_acceptors[acceptor[0]]:
            h_bonds.append(eval_H_bonding.evaluate(donor, acceptor[1], eq_class[1][0], prot_library))
            # if the H-bond evaluation is not successful, print the error message and exit
            successful_completion = h_bonds[-1][0]
            if not successful_completion:
                for note in h_bonds[-1][1]:
                    print(note)
                sys.exit(1)
        acceptor_h_bonds.append(h_bonds)
    # extract the atoms near the protonated donors of interest that can act as H-bond acceptors
    acceptors_near_prot_donors = []
    for atom_group in enumerate(atoms_near_prot_donors):
        list_of_acceptors = []
        for atom in atom_group[1]:
            # the side chain acceptor atoms in ASN, GLN, and HIS residues will be collected separately
            if not (atom[1] != "O" and atom[2] in ["ASN", "GLN", "HIS"]):
                # the acceptors should not belong to the nucleobase of the donor
                if not ((prot_donor_list[atom_group[0]][2] == atom[2] == "A" and
                         prot_donor_list[atom_group[0]][3] == atom[3] and
                         prot_donor_list[atom_group[0]][4] == atom[4] and
                         atom[1] in adenine_nuc_acceptors) or
                        (prot_donor_list[atom_group[0]][2] == atom[2] == "C" and
                         prot_donor_list[atom_group[0]][3] == atom[3] and
                         prot_donor_list[atom_group[0]][4] == atom[4] and
                         atom[1] in cytosine_nuc_acceptors) or
                        (prot_donor_list[atom_group[0]][2] == atom[2] == "G" and
                         prot_donor_list[atom_group[0]][3] == atom[3] and
                         prot_donor_list[atom_group[0]][4] == atom[4] and
                         atom[1] in guanine_nuc_acceptors) or
                        (prot_donor_list[atom_group[0]][2] == atom[2] == "U" and
                         prot_donor_list[atom_group[0]][3] == atom[3] and
                         prot_donor_list[atom_group[0]][4] == atom[4] and
                         atom[1] in uracil_nuc_acceptors)):
                    for acceptor in acceptor_atoms:
                        if atom[1] == acceptor[1] and atom[2] == acceptor[0]:
                            list_of_acceptors.append(atom)
        acceptors_near_prot_donors.append(list_of_acceptors)
    # acquire the H-bonding geometry measurements for all acceptors near each protonated donor
    prot_donor_h_bonds = []
    for donor in enumerate(prot_donor_list):
        h_bonds = []
        for acceptor in acceptors_near_prot_donors[donor[0]]:
            h_bonds.append(eval_H_bonding.evaluate(donor[1], acceptor, eq_class[1][0], prot_library))
            # if the H-bond evaluation is not successful, print the error message and exit
            successful_completion = h_bonds[-1][0]
            if not successful_completion:
                for note in h_bonds[-1][1]:
                    print(note)
                sys.exit(1)
        prot_donor_h_bonds.append(h_bonds)
    # using a greater search distance for the side chains of ASN, GLN, and HIS residues, store a list of atoms near the
    # acceptors of interest
    # additionally, construct a list of atoms near the protonated donors of interest using the greater search distance
    atoms_near_acceptors_amb = []
    atoms_near_prot_donors_amb = []
    for acceptor in stored.acceptor_list:
        stored.nearby_atoms = []
        cmd.iterate(f'index {acceptor[0]} around {search_dist_amb}',
                    'stored.nearby_atoms.append((index, name, resn, resi, chain))')
        atoms_near_acceptors_amb.append(stored.nearby_atoms)
        if (acceptor[2], acceptor[1]) in prot_donors_of_interest:
            atoms_near_prot_donors_amb.append(stored.nearby_atoms)
    # extract the atoms that can act as H-bond donors that belong to the side chains of ASN, GLN, and HIS residues
    donors_near_acceptors_amb = []
    for atom_group in enumerate(atoms_near_acceptors_amb):
        list_of_donors = []
        for atom in atom_group[1]:
            # only extract side chain donor atoms in ASN, GLN, and HIS residues
            if atom[1] != "N" and atom[2] in ["ASN", "GLN", "HIS"]:
                for donor in donor_atoms:
                    if atom[1] == donor[1] and atom[2] == donor[0]:
                        list_of_donors.append(atom)
        donors_near_acceptors_amb.append(list_of_donors)
    # acquire the H-bonding geometry measurements for side chain donor atoms in ASN, GLN, and HIS residues near each
    # acceptor
    acceptor_h_bonds_amb = []
    for acceptor in enumerate(stored.acceptor_list):
        h_bonds = []
        for donor in donors_near_acceptors_amb[acceptor[0]]:
            h_bonds.append(eval_H_bonding.evaluate(donor, acceptor[1], eq_class[1][0], prot_library))
            # if the H-bond evaluation is not successful, print the error message and exit
            successful_completion = h_bonds[-1][0]
            if not successful_completion:
                for note in h_bonds[-1][1]:
                    print(note)
                sys.exit(1)
        acceptor_h_bonds_amb.append(h_bonds)
    # extract the atoms near the protonated donors of interest that can act as H-bond acceptors that belong to the
    # side chains of ASN, GLN, and HIS residues
    acceptors_near_prot_donors_amb = []
    for atom_group in enumerate(atoms_near_prot_donors_amb):
        list_of_acceptors = []
        for atom in atom_group[1]:
            # only extract side chain acceptor atoms in ASN, GLN, and HIS residues
            if atom[1] != "O" and atom[2] in ["ASN", "GLN", "HIS"]:
                for acceptor in acceptor_atoms:
                    if atom[1] == acceptor[1] and atom[2] == acceptor[0]:
                        list_of_acceptors.append(atom)
        acceptors_near_prot_donors_amb.append(list_of_acceptors)
    # acquire the H-bonding geometry measurements for side chain acceptor atoms in ASN, GLN, and HIS residues near each
    # protonated donor
    prot_donor_h_bonds_amb = []
    for donor in enumerate(prot_donor_list):
        h_bonds = []
        for acceptor in acceptors_near_prot_donors_amb[donor[0]]:
            h_bonds.append(eval_H_bonding.evaluate(donor[1], acceptor, eq_class[1][0], prot_library))
            # if the H-bond evaluation is not successful, print the error message and exit
            successful_completion = h_bonds[-1][0]
            if not successful_completion:
                for note in h_bonds[-1][1]:
                    print(note)
                sys.exit(1)
        prot_donor_h_bonds_amb.append(h_bonds)
    # construct a list of all nucleobases containing all donors and acceptors of interest
    nucleobase_list = []
    for i in range(len(stored.donor_list)):
        if (stored.donor_list[i][2] == stored.acceptor_list[i*2][2] == stored.acceptor_list[i*2+1][2] and
                stored.donor_list[i][3] == stored.acceptor_list[i*2][3] == stored.acceptor_list[i*2+1][3] and
                stored.donor_list[i][4] == stored.acceptor_list[i*2][4] == stored.acceptor_list[i*2+1][4]):
            nucleobase_list.append((stored.donor_list[i][2]+stored.donor_list[i][3]+stored.donor_list[i][4],
                                    stored.donor_list[i][0], stored.donor_list[i][1],
                                    stored.acceptor_list[i*2][0], stored.acceptor_list[i*2][1],
                                    stored.acceptor_list[i*2+1][0], stored.acceptor_list[i*2+1][1]))
        else:
            print("Error: The code was unable to complete construction of nucleobase_list.")
            sys.exit(1)

    # TODO collect and report b-factors

    with open(f"{eq_class[0]}_nuc_test_data.csv", "w") as csv_file:
        writer = csv.writer(csv_file)
        for nuc in nucleobase_list:
            writer.writerow(nuc)

    with open(f"{eq_class[0]}_hbond_test_data.csv", "w") as csv_file:
        writer = csv.writer(csv_file)
        for x in enumerate(donor_h_bonds):
            for y in enumerate(x[1]):
                for z in y[1][1]:
                    row = []
                    row.extend(stored.donor_list[x[0]])
                    row.extend(acceptors_near_donors[x[0]][y[0]])
                    row.extend(z)
                    writer.writerow(row)
        for x in enumerate(donor_h_bonds_amb):
            for y in enumerate(x[1]):
                for z in y[1][1]:
                    row = []
                    row.extend(stored.donor_list[x[0]])
                    row.extend(acceptors_near_donors_amb[x[0]][y[0]])
                    row.extend(z)
                    writer.writerow(row)
        for x in enumerate(acceptor_h_bonds):
            for y in enumerate(x[1]):
                for z in y[1][1]:
                    row = []
                    row.extend(donors_near_acceptors[x[0]][y[0]])
                    row.extend(stored.acceptor_list[x[0]])
                    row.extend(z)
                    writer.writerow(row)
        for x in enumerate(acceptor_h_bonds_amb):
            for y in enumerate(x[1]):
                for z in y[1][1]:
                    row = []
                    row.extend(donors_near_acceptors_amb[x[0]][y[0]])
                    row.extend(stored.acceptor_list[x[0]])
                    row.extend(z)
                    writer.writerow(row)
        for x in enumerate(prot_donor_h_bonds):
            for y in enumerate(x[1]):
                for z in y[1][1]:
                    row = []
                    row.extend(prot_donor_list[x[0]])
                    row.extend(acceptors_near_prot_donors[x[0]][y[0]])
                    row.extend(z)
                    writer.writerow(row)
        for x in enumerate(prot_donor_h_bonds_amb):
            for y in enumerate(x[1]):
                for z in y[1][1]:
                    row = []
                    row.extend(prot_donor_list[x[0]])
                    row.extend(acceptors_near_prot_donors_amb[x[0]][y[0]])
                    row.extend(z)
                    writer.writerow(row)

end = datetime.now()
print(end - start)

#     cmd.save(f'{modified_mmCIF_directory}/{eq_class[0]}.cif')
#     cmd.delete('all')
