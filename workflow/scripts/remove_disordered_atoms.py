"""
This module defines a function named remove that uses the Biopython library to load an mmCIF file and then save the
loaded structure with only one set of the disordered atoms, resulting in a "no disorder" structure that only contains a
single set of atoms for each residue.
"""

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO, Select
from Bio.PDB.Selection import unfold_entities


def remove(dir_of_disordered_mmcif, eq_class_mem_id):

    # Define a subclass of Bio.PDB.PDBIO.Select to select the atoms.
    class AtomSelect(Select):
        def __init__(self, list_of_atoms):
            self.list_of_atoms = list_of_atoms
        def accept_atom(self, atm):
            if atm.get_full_id() in self.list_of_atoms:
                self.list_of_atoms.remove(atm.get_full_id())
                return True
            else:
                return False

    # Create a Biopython structure object.
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(eq_class_mem_id, dir_of_disordered_mmcif + eq_class_mem_id + "_with_disorder.cif")

    # Create a list of atoms to keep in the output structure.
    atom_list = []
    for atom in unfold_entities(structure, 'A'):
        atom_list.append(atom.get_full_id())

    # Write the fragments to CIF files.
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(dir_of_disordered_mmcif + eq_class_mem_id + "_no_disorder.cif", AtomSelect(atom_list))
