"""
The sasa function provided by this module takes a list of amine donors and calculates the solvent accessible surface
area (SASA) of each amine. If the encompassing structure contains hydrogens, it is saved to a new CIF file without the
hydrogens, the new structure is loaded, then the SASA calculations are performed.
"""

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.mmcifio import MMCIFIO, Select
from Bio.PDB import Selection


def sasa(pdb_id, model, eq_class_mem_id, donor_list, original_mmcif_dir, h_omitted_mmcif_dir):

    # Initialize the mmCIF parser and load the original (unmodified) mmCIF file.
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(pdb_id, f"{original_mmcif_dir}{pdb_id.lower()}.cif")

    # Check if there are hydrogens in the structure.
    h_present = False
    for atom in Selection.unfold_entities(structure, "A"):
        if atom.element == "H":
            h_present = True
            break

    # Define a function to calculate the SASA values of the amines.
    def calculate(donors, bio_pdb_structure, pymol_model):
        sr = ShrakeRupley()
        sr.compute(bio_pdb_structure)
        sasa_list = []
        for pymol_amine in donors:
            chain = pymol_amine[4]
            hetfield = " "
            # Residue numbers and insertion codes are reported together in pymol_amine[3]. If an insertion code is present,
            # tease it apart from the residue number.
            if pymol_amine[3][-1].isdigit():
                resseq = int(pymol_amine[3])
                icode = " "
            else:
                # This module does not currently handle insertion codes that are more than one character in length.
                if not pymol_amine[3][-2].isdigit():
                    return ("Error: There is an insertion code that is more than one character in length in "
                            f"{eq_class_mem_id}.")
                resseq = int(pymol_amine[3][:-1])
                icode = pymol_amine[3][-1]
            name = pymol_amine[1]
            # PyMOL represents a blank altloc as "" while Bio.PDB represents a blank altloc as " ".
            if pymol_amine[6] == "":
                altloc = " "
            else:
                altloc = pymol_amine[6]
            bio_pdb_amine = bio_pdb_structure[int(pymol_model)-1][chain][(hetfield, resseq, icode)][name]
            # If the alternative location of bio_pdb_amine matches that of pymol_amine, get the SASA of the atom. If
            # not, select the correct alternative location and then get the SASA.
            if bio_pdb_amine.get_altloc() == altloc:
                sasa_list.append(bio_pdb_amine.sasa)
            else:
                bio_pdb_amine.disordered_select(altloc)
                sasa_list.append(bio_pdb_amine.sasa)
        return sasa_list

    # If there are no hydrogens, proceed with the SASA calculations.
    if not h_present:
        sasa_values = calculate(donor_list, structure, model)
        return sasa_values

    # Otherwise, save the structure without hydrogens, load the new structure, and then perform the SASA calculations.
    else:
        # Define a subclass of Bio.PDB.PDBIO.Select to select non-hydrogen atoms.
        class AtomSelect(Select):
            def accept_atom(self, atom):
                if atom.element == "H":
                    return False
                else:
                    return True

        # Write the structures to mmCIF files with hydrogens omitted.
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(f"{h_omitted_mmcif_dir}{eq_class_mem_id}.cif", AtomSelect())

        # Load the mmCIF file with hydrogens omitted.
        new_structure = parser.get_structure(f"h_omitted_{pdb_id}", f"{h_omitted_mmcif_dir}{eq_class_mem_id}.cif")

        # Now, proceed with the SASA calculations.
        sasa_values = calculate(donor_list, new_structure, model)
        return sasa_values
