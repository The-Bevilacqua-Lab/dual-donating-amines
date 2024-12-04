"""
Create the atom_depth function, which uses the Bio.PDB library to calculate the minimum distance between each atom
within a list of donors and the surface of the structure.
"""
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.ResidueDepth import get_surface, min_dist


def atom_depth(pdb_id, model, eq_class_mem_id, donor_list):
    parser = MMCIFParser()
    structure = parser.get_structure(pdb_id, f"{snakemake.config['original_mmcif_dir']}{pdb_id.lower()}.cif")
    surface = get_surface(structure[model])
    distances = []
    for pymol_donor in donor_list:
        chain = pymol_donor[4]
        hetfield = " "
        # Residue numbers and insertion codes are reported together in pymol_donor[3]. If an insertion code is present,
        # tease it apart from the residue number.
        if pymol_donor[3][-1].isdigit():
            resseq = pymol_donor[3]
            icode = " "
        else:
            resseq = pymol_donor[:-1]
            icode = pymol_donor[-1]
            # This module does not currently handle insertion codes that are more than one character in length.
            if pymol_donor[3][-2].isdigit():
                return ("Error: There is an insertion code that is more than one character in length in "
                        f"{eq_class_mem_id}.")
        name = pymol_donor[1]
        # PyMOL represents a blank altloc as "" while Bio.PDB represents a blank altloc as " ".
        if pymol_donor[6] == "":
            altloc = " "
        else:
            altloc = pymol_donor[6]
        bio_pdb_donor = structure[model][chain][(hetfield, resseq, icode)][name]
        # If the alternative location of bio_pdb_donor matches that of pymol_donor, get the coordinates of the atom. If
        # not, select the correct alternative location and then get the coordinates.
        if bio_pdb_donor.get_altloc() == altloc:
            coord = bio_pdb_donor.get_coord()
        else:
            bio_pdb_donor.disordered_select(altloc)
            coord = bio_pdb_donor.get_coord()
        distances.append(min_dist(coord, surface))
    return distances
