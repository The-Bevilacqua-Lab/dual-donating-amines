from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley


def sasa(pdb_id, model, eq_class_mem_id, donor_list, original_mmcif_dir):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(pdb_id, f"{original_mmcif_dir}{pdb_id.lower()}.cif")
    sr = ShrakeRupley()
    sr.compute(structure)
    sasa_list = []
    for pymol_amine in donor_list:
        chain = pymol_amine[4]
        hetfield = " "
        # Residue numbers and insertion codes are reported together in pymol_amine[3]. If an insertion code is present,
        # tease it apart from the residue number.
        if pymol_amine[3][-1].isdigit():
            resseq = int(pymol_amine[3])
            icode = " "
        else:
            resseq = int(pymol_amine[3][:-1])
            icode = pymol_amine[3][-1]
            # This module does not currently handle insertion codes that are more than one character in length.
            if not pymol_amine[3][-2].isdigit():
                return ("Error: There is an insertion code that is more than one character in length in "
                        f"{eq_class_mem_id}.")
        name = pymol_amine[1]
        # PyMOL represents a blank altloc as "" while Bio.PDB represents a blank altloc as " ".
        if pymol_amine[6] == "":
            altloc = " "
        else:
            altloc = pymol_amine[6]
        bio_pdb_amine = structure[int(model)-1][chain][(hetfield, resseq, icode)][name]
        # If the alternative location of bio_pdb_amine matches that of pymol_amine, get the coordinates of the atom. If
        # not, select the correct alternative location and then get the coordinates.
        if bio_pdb_amine.get_altloc() == altloc:
            sasa_list.append(bio_pdb_amine.sasa)
        else:
            bio_pdb_amine.disordered_select(altloc)
            sasa_list.append(bio_pdb_amine.sasa)
    return sasa_list
