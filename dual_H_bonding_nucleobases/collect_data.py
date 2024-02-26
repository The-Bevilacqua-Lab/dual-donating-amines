"""
This script reads information from the equivalence class file created by parse_nrlist.py. The name of the equivalence
class file must be provided as the first argument, and this file must exist within a folder named eq_class_files in the
current working directory. This script works with PyMOL to save an mmCIF file of the original structure to the
original_mmCIF_files folder and an mmCIF file of the structure with alternative conformations removed and hydrogens
added to the modified_mmCIF_files folder. It works with the modified structure to collect data on potential hydrogen
bonds involving atoms of interest and nucleobase atom b-factors in the representative RNA chains. This data is then
written to three separate csv files that are saved in three separate folders. If commit_hash is provided as an argument,
the commit hash of the repo will also be written within a commented line to the data files if no uncommitted changes
have been made to the repo.
"""

import sys
import os
import csv
import subprocess
from pymol import cmd
from pymol import stored
import numpy as np
from datetime import datetime
import const
import eval_H_bonding
import remove_alt_conf

start = datetime.now()

# check to make sure the folder named eq_class_files and the specified equivalence class file exists
eq_class_file = ""
try:
    if not os.path.isdir(os.getcwd() + "/eq_class_files"):
        print(f"Error: The folder named eq_class_files was not found in the current working directory.")
        sys.exit(1)
    elif not os.path.isfile(os.getcwd() + "/eq_class_files/" + sys.argv[1]):
        print(f"Error: The file named {sys.argv[1]} was not found in the eq_class_files folder.")
        sys.exit(1)
    else:
        eq_class_file = os.getcwd() + "/eq_class_files/" + sys.argv[1]
except IndexError:
    print("Error: The name of the equivalence class file must be provided as the first argument.")
    sys.exit(1)

# iterate through the lines of the equivalence class file and collect the equivalence class name, PDB ID, model info,
# and chain info for the representative structure
with open("eq_class_files/" + eq_class_file, mode='r') as read_file:
    eq_class = []
    for line in csv.reader(read_file):
        if line[0][0] != "#":
            eq_class.append(line)

# check whether an original_mmCIF_files folder exists, and if it does not exist, create the folder
original_mmCIF_dir = os.getcwd() + "/original_mmCIF_files"
if not os.path.isdir(original_mmCIF_dir):
    os.mkdir(original_mmCIF_dir)

# change fetch_path to the original_mmCIF_files folder in the current working directory so that PyMOL drops the mmCIF
# files into this folder when running fetch
cmd.set('fetch_path', cmd.exp_path(original_mmCIF_dir))

# construct two lists of tuples containing atom and residue names that describe either donor or acceptor atoms
# additionally, construct a list of tuples containing atom and residue names that describe rotatable donor atoms
donor_atoms = []
rotatable_donor_atoms = []
acceptor_atoms = []
for residue in const.RESIDUE_LIBRARY:
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
for residue in const.RESIDUE_LIBRARY:
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

# construct four strings that can be used with PyMOL to select all possible donor, all possible protonated donor
# (formal charge of +1), all possible acceptor, or all possible deprotonated acceptor (formal charge of -1) atoms of
# particular interest
donors_of_interest_str = ''
last_atom = const.DONORS_OF_INTEREST[-1]
for atom in const.DONORS_OF_INTEREST:
    donors_of_interest_str += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        donors_of_interest_str += ' '
prot_donors_of_interest_str = ''
last_atom = const.PROT_DONORS_OF_INTEREST[-1]
for atom in const.PROT_DONORS_OF_INTEREST:
    prot_donors_of_interest_str += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        prot_donors_of_interest_str += ' '
acceptors_of_interest_str = ''
last_atom = const.ACCEPTORS_OF_INTEREST[-1]
for atom in const.ACCEPTORS_OF_INTEREST:
    acceptors_of_interest_str += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        acceptors_of_interest_str += ' '
deprot_acceptors_of_interest_str = ''
last_atom = const.DEPROT_ACCEPTORS_OF_INTEREST[-1]
for atom in const.DEPROT_ACCEPTORS_OF_INTEREST:
    deprot_acceptors_of_interest_str += f'(resn {atom[0]} and name {atom[1]})'
    if atom != last_atom:
        deprot_acceptors_of_interest_str += ' '

# define the donor and acceptor atoms of the RNA nucleobases
adenine_nuc_donors = ["N6"]
adenine_nuc_acceptors = ["N1", "N3", "N7"]
cytosine_nuc_donors = ["N4"]
cytosine_nuc_acceptors = ["O2", "N3"]
guanine_nuc_donors = ["N1", "N2"]
guanine_nuc_acceptors = ["N3", "O6", "N7"]
uracil_nuc_donors = ["N3"]
uracil_nuc_acceptors = ["O2", "O4"]

# using the library provided by the const module, create a new library with protonated A, C, and G residues
expanded_library = const.RESIDUE_LIBRARY
for res in expanded_library:
    if res['res'] == 'A':
        res['don'].append(["N1", "N", False, 1, "C2"])
        res['don'].append(["N3", "N", False, 1, "C2"])
        res['don'].append(["N7", "N", False, 1, "C5"])
    if res['res'] == 'C':
        res['don'].append(["N3", "N", False, 1, "C2"])
    if res['res'] == 'G':
        res['don'].append(["N3", "N", False, 1, "C2"])
        res['don'].append(["N7", "N", False, 1, "C5"])
        res['acc'].append(["N1", "N", False, -1, "C2"])
    if res['res'] == 'U':
        res['acc'].append(["N3", "N", False, -1, "C2"])

# define the distance used to search for nearby donor or acceptor atoms
search_dist = 4.1
# define the distance used to search for nearby donor or acceptor atoms that belong to side chains that have
# conformations that may have been ambiguously assigned
search_dist_amb = search_dist + 3.0

# identify the chains of all the representative RNAs for the equivalence class
rep_rna_list = []
for i in range(len(eq_class[2])):
    # the create function below assumes that the different representative RNA chains belong to the same model
    # if this is not the case, print an error message and exit
    if eq_class[2][i] != eq_class[2][0]:
        print(f"Error: There are representative RNA chains that belong to different models for equivalence class "
              f"{eq_class[0][0]}.")
        sys.exit(1)
    rep_rna_list.append(f"chain {eq_class[3][i]}")
rep_rna = " ".join(rep_rna_list)

# retrieve the structure that contains the representative RNA chains
cmd.fetch(eq_class[1][0])

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

# Change the formal charge on ARG(NH2), TYR(OH), HIS(ND1), A(N1), A(N3), C(N3), and G(N3) to +1. With a formal charge
# of +1, PyMOL will add hydrogens to these atoms when the following h_add command is used. The NH2 atoms of ARG residues
# typically have a formal charge of 0. With a formal charge of 0, PyMOL only adds a single hydrogen to these atoms when
# the h_add command is issued using a script. When the h_add command is issued through the PyMOL GUI, PyMOL correctly
# adds two hydrogens without needing to increase the formal charge of the NH2 atoms to +1. This must be some sort of bug
# with PyMOL 2.5.0 Open-Source. Other versions of PyMOL may not have this same behavior.
cmd.alter(f'(resn ARG and name NH2) (resn TYR and name OH) (resn HIS and name ND1) {prot_donors_of_interest_str}',
          'formal_charge=1')

# add hydrogens to non-rotatable donors that are a part of or near the representative RNA
cmd.h_add(f'(({donor_string} {prot_donors_of_interest_str}) and not ({rotatable_donor_string})) within '
          f'{search_dist_amb + 3.0} of ({rep_rna})')

# randomly sample 100 atom indices in the structure and record the info of the associated atoms for later comparison
num_atoms = cmd.count_atoms('all')
indices = np.random.default_rng().integers(low=1, high=num_atoms+1, size=100)
stored.check_one = []
for index in indices:
    cmd.iterate(f'index {index}', 'stored.check_one.append((index, name, resn, resi, chain))')

# store a list of donors of interest from the representative RNA
stored.donor_list = []
cmd.iterate(f'({rep_rna}) and ({donors_of_interest_str})',
            'stored.donor_list.append((index, name, resn, resi, chain))')

# store a list of atoms near the donors of interest
atoms_near_donors = []
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
        h_bonds.append(eval_H_bonding.evaluate(donor[1], acceptor, eq_class[1][0], expanded_library))
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
        h_bonds.append(eval_H_bonding.evaluate(donor[1], acceptor, eq_class[1][0], expanded_library))
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
for acceptor in stored.acceptor_list:
    stored.nearby_atoms = []
    cmd.iterate(f'index {acceptor[0]} around {search_dist}',
                'stored.nearby_atoms.append((index, name, resn, resi, chain))')
    atoms_near_acceptors.append(stored.nearby_atoms)
    if (acceptor[2], acceptor[1]) in const.PROT_DONORS_OF_INTEREST:
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
                if atom[1] == "O3'" or atom[1] == "O5'":
                    if eval_H_bonding.terminal_donor(atom):
                        list_of_donors.append(atom)
                else:
                    for donor in donor_atoms:
                        if atom[1] == donor[1] and atom[2] == donor[0]:
                            list_of_donors.append(atom)
    donors_near_acceptors.append(list_of_donors)

# acquire the H-bonding geometry measurements for all donors near each acceptor
acceptor_h_bonds = []
for acceptor in enumerate(stored.acceptor_list):
    h_bonds = []
    for donor in donors_near_acceptors[acceptor[0]]:
        h_bonds.append(eval_H_bonding.evaluate(donor, acceptor[1], eq_class[1][0], expanded_library))
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
        h_bonds.append(eval_H_bonding.evaluate(donor[1], acceptor, eq_class[1][0], expanded_library))
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
    if (acceptor[2], acceptor[1]) in const.PROT_DONORS_OF_INTEREST:
        atoms_near_prot_donors_amb.append(stored.nearby_atoms)

# extract the atoms that can act as H-bond donors that belong to the side chains of ASN, GLN, and HIS residues
donors_near_acceptors_amb = []
for atom_group in enumerate(atoms_near_acceptors_amb):
    list_of_donors = []
    for atom in atom_group[1]:
        # only extract side chain donor atoms in ASN, GLN, and HIS residues
        if atom[1] != "N" and atom[2] in ["ASN", "GLN", "HIS"]:
            if atom[1] == "O3'" or atom[1] == "O5'":
                if eval_H_bonding.terminal_donor(atom):
                    list_of_donors.append(atom)
            else:
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
        h_bonds.append(eval_H_bonding.evaluate(donor, acceptor[1], eq_class[1][0], expanded_library))
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
        h_bonds.append(eval_H_bonding.evaluate(donor[1], acceptor, eq_class[1][0], expanded_library))
        # if the H-bond evaluation is not successful, print the error message and exit
        successful_completion = h_bonds[-1][0]
        if not successful_completion:
            for note in h_bonds[-1][1]:
                print(note)
            sys.exit(1)
    prot_donor_h_bonds_amb.append(h_bonds)

# store a list of deprotonated acceptors of interest (formal charge of -1) from the representative RNA
stored.deprot_acceptor_list = []
cmd.iterate(f'({rep_rna}) and ({deprot_acceptors_of_interest_str})',
            'stored.deprot_acceptor_list.append((index, name, resn, resi, chain))')

# store a list of atoms near the deprotonated acceptors of interest
atoms_near_deprot_acceptors = []
for acceptor in stored.deprot_acceptor_list:
    stored.nearby_atoms = []
    cmd.iterate(f'index {acceptor[0]} around {search_dist}',
                'stored.nearby_atoms.append((index, name, resn, resi, chain))')
    atoms_near_deprot_acceptors.append(stored.nearby_atoms)

# extract the atoms near the deprotonated acceptors of interest that can act as H-bond donors
donors_near_deprot_acceptors = []
for atom_group in enumerate(atoms_near_deprot_acceptors):
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
                if atom[1] == "O3'" or atom[1] == "O5'":
                    if eval_H_bonding.terminal_donor(atom):
                        list_of_donors.append(atom)
                else:
                    for donor in donor_atoms:
                        if atom[1] == donor[1] and atom[2] == donor[0]:
                            list_of_donors.append(atom)
    donors_near_deprot_acceptors.append(list_of_donors)

# acquire the H-bonding geometry measurements for all donors near each deprotonated acceptor
deprot_acceptor_h_bonds = []
for acceptor in enumerate(stored.deprot_acceptor_list):
    h_bonds = []
    for donor in donors_near_deprot_acceptors[acceptor[0]]:
        h_bonds.append(eval_H_bonding.evaluate(donor, acceptor[1], eq_class[1][0], expanded_library))
        # if the H-bond evaluation is not successful, print the error message and exit
        successful_completion = h_bonds[-1][0]
        if not successful_completion:
            for note in h_bonds[-1][1]:
                print(note)
            sys.exit(1)
    deprot_acceptor_h_bonds.append(h_bonds)

# using a greater search distance for the side chains of ASN, GLN, and HIS residues, store a list of atoms near the
# deprotonated acceptors of interest (formal charge of -1)
atoms_near_deprot_acceptors_amb = []
for acceptor in stored.deprot_acceptor_list:
    stored.nearby_atoms = []
    cmd.iterate(f'index {acceptor[0]} around {search_dist_amb}',
                'stored.nearby_atoms.append((index, name, resn, resi, chain))')
    atoms_near_deprot_acceptors_amb.append(stored.nearby_atoms)

# extract the atoms that can act as H-bond donors that belong to the side chains of ASN, GLN, and HIS residues
donors_near_deprot_acceptors_amb = []
for atom_group in enumerate(atoms_near_deprot_acceptors_amb):
    list_of_donors = []
    for atom in atom_group[1]:
        # only extract side chain donor atoms in ASN, GLN, and HIS residues
        if atom[1] != "N" and atom[2] in ["ASN", "GLN", "HIS"]:
            if atom[1] == "O3'" or atom[1] == "O5'":
                if eval_H_bonding.terminal_donor(atom):
                    list_of_donors.append(atom)
            else:
                for donor in donor_atoms:
                    if atom[1] == donor[1] and atom[2] == donor[0]:
                        list_of_donors.append(atom)
    donors_near_deprot_acceptors_amb.append(list_of_donors)

# acquire the H-bonding geometry measurements for side chain donor atoms in ASN, GLN, and HIS residues near each
# deprotonated acceptor
deprot_acceptor_h_bonds_amb = []
for acceptor in enumerate(stored.deprot_acceptor_list):
    h_bonds = []
    for donor in donors_near_deprot_acceptors_amb[acceptor[0]]:
        h_bonds.append(eval_H_bonding.evaluate(donor, acceptor[1], eq_class[1][0], expanded_library))
        # if the H-bond evaluation is not successful, print the error message and exit
        successful_completion = h_bonds[-1][0]
        if not successful_completion:
            for note in h_bonds[-1][1]:
                print(note)
            sys.exit(1)
    deprot_acceptor_h_bonds_amb.append(h_bonds)

# construct a list of all nucleobases containing donors, protonated donors, acceptors, and deprotonated acceptors of
# interest
nucleobase_list = []
nucleobase_list.extend(stored.donor_list)
nucleobase_list.extend(stored.acceptor_list)
nucleobase_list.extend(stored.deprot_acceptor_list)

# collect the b-factors of nucleobase atoms within the representative structure
stored.res_list = []
cmd.iterate(f"name C1' and ({rep_rna})", "stored.res_list.append((index, name, resn, resi, chain))")
nuc_b_factors = []
for res in stored.res_list:
    stored.b_factors = []
    cmd.iterate(f"(byres index {res[0]}) and sidechain", "stored.b_factors.append(b)")
    nuc_b_factors.append([f"{res[2]}{res[3]}{res[4]}"])
    nuc_b_factors[-1].extend(stored.b_factors)

# record info on the atoms associated with the previously determined 100 atom indices and check whether there are any
# changes
stored.check_two = []
for index in indices:
    cmd.iterate(f'index {index}', 'stored.check_two.append((index, name, resn, resi, chain))')
if not stored.check_one == stored.check_two:
    print(f"Error: The indices in PDB ID {eq_class[1][0]} changed.")
    sys.exit(1)

# If commit_hash is supplied as the first argument, check if any changes have been made to the repo and get the hash
# of the current git commit. If uncommitted changes have been made, print an error message and exit.
repo_changes = ""
commit_hash = ""
if "commit_hash" in sys.argv:
    repo_changes = subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=no"],
                                           cwd=os.path.dirname(os.path.realpath(__file__))).decode('ascii').strip()
    if not repo_changes:
        commit_hash = subprocess.check_output(["git", "rev-parse", "HEAD"],
                                              cwd=os.path.dirname(os.path.realpath(__file__))).decode('ascii').strip()
    else:
        print(f"Error: Uncommitted changes have been made to the repo.")
        sys.exit(1)

# check whether a b_factor_data folder exists, and if it does not exist, create the folder
b_factor_data_dir = os.getcwd() + "/b_factor_data"
if not os.path.isdir(b_factor_data_dir):
    os.mkdir(b_factor_data_dir)

# write a csv containing the b-factors of all nucleobase atoms in the representative structure
if not os.path.isfile(f"{b_factor_data_dir}/{eq_class[0][0]}_b_factor_data.csv"):
    with open(f"{b_factor_data_dir}/{eq_class[0][0]}_b_factor_data.csv", "w") as csv_file:
        writer = csv.writer(csv_file)
        if commit_hash:
            writer.writerow([f"# dual-H-bonding-nucleobases repo git commit hash: {commit_hash}"])
        writer.writerow([f"# input file: {eq_class_file}"])
        writer.writerow([f"# file created on: {datetime.now().strftime('%y-%m-%d %H:%M:%S.%f')}"])
        for nuc in nuc_b_factors:
            writer.writerow(nuc)
else:
    print(f"Error: The file named {eq_class[0][0]}_b_factor_data.csv already exists.")
    sys.exit(1)

# check whether a nuc_data folder exists, and if it does not exist, create the folder
nuc_data_dir = os.getcwd() + "/nuc_data"
if not os.path.isdir(nuc_data_dir):
    os.mkdir(nuc_data_dir)

# write a csv containing all nucleobases that contain donors, protonated donors, acceptors, and deprotonated acceptors
# of interest
if not os.path.isfile(f"{nuc_data_dir}/{eq_class[0][0]}_nuc_data.csv"):
    with open(f"{nuc_data_dir}/{eq_class[0][0]}_nuc_data.csv", "w") as csv_file:
        writer = csv.writer(csv_file)
        if commit_hash:
            writer.writerow([f"# dual-H-bonding-nucleobases repo git commit hash: {commit_hash}"])
        writer.writerow([f"# input file: {eq_class_file}"])
        writer.writerow([f"# file created on: {datetime.now().strftime('%y-%m-%d %H:%M:%S.%f')}"])
        for nuc in nucleobase_list:
            writer.writerow(nuc)
else:
    print(f"Error: The file named {eq_class[0][0]}_nuc_data.csv already exists.")
    sys.exit(1)

# check whether a hbond_data folder exists, and if it does not exist, create the folder
hbond_data_dir = os.getcwd() + "/hbond_data"
if not os.path.isdir(hbond_data_dir):
    os.mkdir(hbond_data_dir)

# write a csv containing all H-bonding information
if not os.path.isfile(f"{hbond_data_dir}/{eq_class[0][0]}_hbond_data.csv"):
    with open(f"{hbond_data_dir}/{eq_class[0][0]}_hbond_data.csv", "w") as csv_file:
        writer = csv.writer(csv_file)
        if commit_hash:
            writer.writerow([f"# dual-H-bonding-nucleobases repo git commit hash: {commit_hash}"])
        writer.writerow([f"# input file: {eq_class_file}"])
        writer.writerow([f"# file created on: {datetime.now().strftime('%y-%m-%d %H:%M:%S.%f')}"])
        for i, don in enumerate(donor_h_bonds):
            for j, acc in enumerate(don):
                for instance in acc[1]:
                    row = []
                    row.extend(stored.donor_list[i])
                    row.extend(acceptors_near_donors[i][j])
                    row.extend(instance)
                    writer.writerow(row)
        for i, don in enumerate(donor_h_bonds_amb):
            for j, acc in enumerate(don):
                for instance in acc[1]:
                    row = []
                    row.extend(stored.donor_list[i])
                    row.extend(acceptors_near_donors_amb[i][j])
                    row.extend(instance)
                    writer.writerow(row)
        for i, acc in enumerate(acceptor_h_bonds):
            for j, don in enumerate(acc):
                for instance in don[1]:
                    row = []
                    row.extend(donors_near_acceptors[i][j])
                    row.extend(stored.acceptor_list[i])
                    row.extend(instance)
                    writer.writerow(row)
        for i, acc in enumerate(acceptor_h_bonds_amb):
            for j, don in enumerate(acc):
                for instance in don[1]:
                    row = []
                    row.extend(donors_near_acceptors_amb[i][j])
                    row.extend(stored.acceptor_list[i])
                    row.extend(instance)
                    writer.writerow(row)
        for i, don in enumerate(prot_donor_h_bonds):
            for j, acc in enumerate(don):
                for instance in acc[1]:
                    row = []
                    row.extend(prot_donor_list[i])
                    row.extend(acceptors_near_prot_donors[i][j])
                    row.extend(instance)
                    writer.writerow(row)
        for i, don in enumerate(prot_donor_h_bonds_amb):
            for j, acc in enumerate(don):
                for instance in acc[1]:
                    row = []
                    row.extend(prot_donor_list[i])
                    row.extend(acceptors_near_prot_donors_amb[i][j])
                    row.extend(instance)
                    writer.writerow(row)
        for i, acc in enumerate(deprot_acceptor_h_bonds):
            for j, don in enumerate(acc):
                for instance in don[1]:
                    row = []
                    row.extend(donors_near_deprot_acceptors[i][j])
                    row.extend(stored.deprot_acceptor_list[i])
                    row.extend(instance)
                    writer.writerow(row)
        for i, acc in enumerate(deprot_acceptor_h_bonds_amb):
            for j, don in enumerate(acc):
                for instance in don[1]:
                    row = []
                    row.extend(donors_near_deprot_acceptors_amb[i][j])
                    row.extend(stored.deprot_acceptor_list[i])
                    row.extend(instance)
                    writer.writerow(row)
else:
    print(f"Error: The file named {eq_class[0][0]}_hbond_data.csv already exists.")
    sys.exit(1)

end = datetime.now()
print(end - start)

# check whether a modified_mmCIF_files folder exists, and if it does not exist, create the folder
modified_mmCIF_dir = os.getcwd() + "/modified_mmCIF_files"
if not os.path.isdir(modified_mmCIF_dir):
    os.mkdir(modified_mmCIF_dir)

cmd.save(f'{modified_mmCIF_dir}/{eq_class[0][0]}.cif')
