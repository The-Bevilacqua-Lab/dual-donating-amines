#################################### Importing libraries ####################################
import os
import json
import requests
import numpy as np
import pandas as pd

from math import floor
from collections import Counter

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

from Bio.PDB import MMCIFParser, Superimposer

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from sklearn.manifold import MDS

import time
from gemmi import cif

#################################### Functions ####################################

# fixing '5E54' becoming '5.00E+54' issue
def fix_PDB_ID(df):
    for i, j in enumerate(df['pdb']):
        if len(str(j)) > 4:
            if '-' in str(j):
                X = ''.join(str(j).split('-')).upper()
                df.loc[i, 'pdb'] = X
            else:
                if '.' in str(j):
                    x = str(j)
                    x1 = x.split('.')[0]
                    x2 = x.split('+')[1]
                    X = x1 + 'E' + x2
                    df.loc[i, 'pdb'] = X
                else:
                    x = str(j)
                    x1 = x.split('+')[0]
                    x2 = x.split('+')[1]
                    X = x1 + x2
                    df.loc[i, 'pdb'] = X
    return df


# functions to generate clipped structures

#function to extract structure file in cif format and convert as a pandas dataframe
def extract_cif(p):

    # Use GEMMI to create a dictionary containing the info in the _atom_site section.
    filename= p.lower()+'.cif'
    mm_cif = cif.read_file(f"resources/original_mmCIF_files/{filename}")
    block = mm_cif.sole_block()

    c_intrst = {'_atom_site.group_PDB': [], '_atom_site.id': [], '_atom_site.type_symbol': [],
                '_atom_site.label_atom_id': [], '_atom_site.label_alt_id': [],
                '_atom_site.label_comp_id': [], '_atom_site.label_asym_id': [],
                '_atom_site.label_entity_id': [], '_atom_site.label_seq_id': [],
                '_atom_site.pdbx_PDB_ins_code': [], '_atom_site.Cartn_x': [],
                '_atom_site.Cartn_y': [], '_atom_site.Cartn_z': [], '_atom_site.occupancy': [],
                '_atom_site.B_iso_or_equiv': [], '_atom_site.pdbx_formal_charge': [],
                '_atom_site.auth_seq_id': [], '_atom_site.auth_comp_id': [],
                '_atom_site.auth_asym_id': [], '_atom_site.auth_atom_id': [],
                '_atom_site.pdbx_PDB_model_num': []}

    table = block.find(list(c_intrst.keys()))
    for row in table:
        new_row = []
        for value in row:
            if "'" in value:
                new_row.append(value.replace("'", "\\'"))
            else:
                new_row.append(value)
        for key, value in zip(c_intrst, new_row):
            c_intrst[key].append(value)

    # Append a row to the end to match what was created by the previous code.
    for key, value in zip(c_intrst, [""] + [None] * 20):
        c_intrst[key].append(value)

    # Create and return a DataFrame based on the dictionary.
    df = pd.DataFrame(c_intrst)
    return df


def clip(p, r, h, F, dr):
    #p= pdb_ID, entries can be case-insensitive
    #r= list of residues 
    #each residue information in the list above will be formatted like: chain_ID.residue_ID
    #h= 'N' or 'Y', if 'N', that means no hydrogen coordinates should be included considering hydrogen coordinates are provided 
    #if 'Y', that means no hydrogen coordinates should be removed if any hydrogen coordinates are provided
    #F= desired output filename (with the extension of '.cif' or '.pdb')
    #dr= path to the clipped structures
    P = p.upper()
    c_ids = list(set([str(i).split('_')[0] for i in r]))
    print("Chains in motif:", c_ids)

    p1 = extract_cif(P)
    if isinstance(p1, str) and p1 == '0':
        return

    if len(p1) > 1:
        p1['temp_col'] = p1[['_atom_site.auth_asym_id', '_atom_site.auth_seq_id', '_atom_site.pdbx_PDB_ins_code']] \
            .fillna("").agg('_'.join, axis=1).str.strip()
        # convert '?' to '*'
        p1['temp_col'] = p1['temp_col'].str.replace('?', '*')
        # instead of filtering out by chain ID, segment ID, residue index, residue ID one after another
        # temp column will have chain ID, residue index, and insertion code all together
        # each row will have their unique items for this column
        p2 = p1[p1['temp_col'].isin(r)]
        p2.index = np.arange(0, len(p2))

        #filtering out or not hydrogen atoms
        if h == 'N':
            p3 = p2[p2['_atom_site.type_symbol'] != 'H']
            p3.index = np.arange(0, len(p3))
        else:
            p3 = p2.copy()

        p3 = p3.drop(columns=['temp_col'])

        #cleaning name for the sugar atoms
        for i, j in enumerate(p3['_atom_site.label_atom_id']):
            if "\\'" in str(j):
                x = str(j).replace("\\'", "'")
                y = x.replace('"', '')
                p3.loc[i, '_atom_site.label_atom_id'] = y
                p3.loc[i, '_atom_site.auth_atom_id'] = y

        if F.split('.')[-1].lower() == 'pdb':
            #converting to pdb
            columns_needed = [
                "_atom_site.id", "_atom_site.label_atom_id", "_atom_site.label_comp_id",
                "_atom_site.auth_asym_id", "_atom_site.auth_seq_id", "_atom_site.pdbx_PDB_ins_code",
                "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z",
                "_atom_site.occupancy", "_atom_site.B_iso_or_equiv", "_atom_site.type_symbol"
            ]
            # Ensure all required columns exist
            for col in columns_needed:
                if col not in p3.columns:
                    p3[col] = ""
            # Open the output PDB file for writing
            with open(dr + F, "w") as pdb_file:
                atom_serial_number = 1
                common_chain_ID = 'A'  #temporary solution
                #main problem is what should be the solution if residues are from multiple chains?
                for _, row in p3.iterrows():
                    # Handle missing insertion codes (CIF uses '?' or '.' for no insertion code)
                    ins_code = row["_atom_site.pdbx_PDB_ins_code"]
                    ins_code = ins_code if ins_code not in ["?", "."] else " "
                    pdb_file.write(
                        f"ATOM  {atom_serial_number:>5}  {row['_atom_site.label_atom_id']:<4}"
                        f"{row['_atom_site.label_comp_id']:>3} {common_chain_ID}"
                        f"{int(row['_atom_site.label_seq_id']):>4}{ins_code:<1}   "
                        f"{float(row['_atom_site.Cartn_x']):>8.3f}{float(row['_atom_site.Cartn_y']):>8.3f}{float(row['_atom_site.Cartn_z']):>8.3f}"
                        f"{float(row['_atom_site.occupancy']):>6.2f}{float(row['_atom_site.B_iso_or_equiv']):>6.2f}          {row['_atom_site.type_symbol']:>2}\n"
                    )
                    atom_serial_number += 1

        else:
            #converting to cif
            #Creating dataframe for the column name
            c_intrst_1 = [
                '_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol',
                '_atom_site.label_atom_id', '_atom_site.label_alt_id',
                '_atom_site.label_comp_id', '_atom_site.label_asym_id',
                '_atom_site.label_entity_id', '_atom_site.label_seq_id',
                '_atom_site.pdbx_PDB_ins_code', '_atom_site.Cartn_x',
                '_atom_site.Cartn_y', '_atom_site.Cartn_z', '_atom_site.occupancy',
                '_atom_site.B_iso_or_equiv', '_atom_site.pdbx_formal_charge',
                '_atom_site.auth_seq_id', '_atom_site.auth_comp_id',
                '_atom_site.auth_asym_id', '_atom_site.auth_atom_id',
                '_atom_site.pdbx_PDB_model_num'
            ]
            c_intrst_2 = [
                'data_', '#', 'loop_',
                '_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol',
                '_atom_site.label_atom_id', '_atom_site.label_alt_id',
                '_atom_site.label_comp_id', '_atom_site.label_asym_id',
                '_atom_site.label_entity_id', '_atom_site.label_seq_id',
                '_atom_site.pdbx_PDB_ins_code', '_atom_site.Cartn_x',
                '_atom_site.Cartn_y', '_atom_site.Cartn_z', '_atom_site.occupancy',
                '_atom_site.B_iso_or_equiv', '_atom_site.pdbx_formal_charge',
                '_atom_site.auth_seq_id', '_atom_site.auth_comp_id',
                '_atom_site.auth_asym_id', '_atom_site.auth_atom_id',
                '_atom_site.pdbx_PDB_model_num'
            ]
            c_intrst_2[0] = 'data_' + P

            p4 = pd.DataFrame(columns=c_intrst_1)
            p4['_atom_site.group_PDB'] = c_intrst_2#this p4 dataframe above p3 will be joined in p4. p3 will help pymol to visualize the file
            p5 = pd.concat([p4, p3], axis=0)
            p5.to_csv(dr + F, header=False, sep=' ', index=False)

# --------- CSV-order-based correspondence for RMSD (NEW CORE) ---------

def _norm_ins_code(x):
    if x is None:
        return ""
    x = str(x).strip()
    if x in ["", " ", ".", "?", "*"]:
        return ""
    return x

def _parse_res_key(res_key):
    #res_key is underscore-joined:
    #  chain_authSeq_insCode
    #e.g. "A_105_*" or "B_47_A" or "A_7_" etc
    parts = str(res_key).split("_")
    while len(parts) < 3:
        parts.append("")
    chain = parts[0].strip()
    resseq = parts[1].strip()
    icode = _norm_ins_code(parts[2])
    return chain, resseq, icode

def build_residue_lookup(structure):
    #Lookup by (chain_id, resseq, icode) based on parsed structure.
    #This is only used within each clipped file to find residues we clipped.
    
    lookup = {}
    for model in structure:
        for chain in model:
            chain_id = str(chain.id).strip()
            for residue in chain:
                hetflag, resseq, icode = residue.id
                icode = _norm_ins_code(icode)
                lookup[(chain_id, str(resseq), icode)] = residue
    return lookup

# will work across sequence differences
BACKBONE_ATOMS = ["P", "C4'"]

PURINE_EXTRA = ["N9", "C6", "C2"]
PYRIM_EXTRA  = ["N1", "C2", "C4"]

def choose_atoms_for_resname(resname, mode="backbone"):
    rn = (resname or "").strip().upper()
    mod_map = {
        "OMG":"G","2MG":"G","G2M":"G",
        "M7A":"A","A2M":"A",
        "5MC":"C","C2M":"C","4AC":"C",
        "OMU":"U","U2M":"U","PSU":"U",
    }
    rn = mod_map.get(rn, rn)

    if mode == "backbone":
        return set(BACKBONE_ATOMS)

    atoms = set(BACKBONE_ATOMS)
    if rn in ["A","G"]:
        atoms |= set(PURINE_EXTRA)
    elif rn in ["C","U"]:
        atoms |= set(PYRIM_EXTRA)
    return atoms

def atom_map_from_csv_order(structure, csv_res_list, mode="backbone"):
    #Build (pos, atom_name) -> Atom map where pos is defined by CSV order.
    lookup = build_residue_lookup(structure)
    amap = {}

    for pos, res_key in enumerate(csv_res_list):
        chain, resseq, icode = _parse_res_key(res_key)
        residue = lookup.get((chain, resseq, icode), None)
        if residue is None:
            continue
        wanted = choose_atoms_for_resname(residue.get_resname(), mode=mode)
        for atom in residue.get_atoms():
            aname = atom.get_name().strip().upper()
            if aname in wanted:
                amap[(pos, aname)] = atom
    return amap

def matched_atom_lists_from_csv(structure1, structure2, csv_list1, csv_list2,
                                mode="backbone", min_common_atoms=6):
    #Matched atom lists using CSV-defined motif order.

    m1 = atom_map_from_csv_order(structure1, csv_list1, mode=mode)
    m2 = atom_map_from_csv_order(structure2, csv_list2, mode=mode)

    common_keys = sorted(set(m1).intersection(m2))
    if len(common_keys) < min_common_atoms:
        raise ValueError(f"Too few common atoms to align: {len(common_keys)}")

    atoms1 = [m1[k] for k in common_keys]
    atoms2 = [m2[k] for k in common_keys]
    return atoms1, atoms2

def align_by_csv_order(cif1, cif2, csv_list1, csv_list2, dr,
                       mode="backbone", min_common_atoms=6):
    parser = MMCIFParser(QUIET=True)
    s1 = parser.get_structure("S1", dr + cif1)
    s2 = parser.get_structure("S2", dr + cif2)

    atoms1, atoms2 = matched_atom_lists_from_csv(
        s1, s2, csv_list1, csv_list2,
        mode=mode, min_common_atoms=min_common_atoms
    )
    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)
    return float(sup.rms)

def align_all_against_all(dr, res_order_map,
                          mode="backbone", min_common_atoms=6, decimals=None):
    #Compute all-vs-all RMSD matrix using CSV-defined residue order.
    #Assumes clipped files are named as <str_ID>.cif

    filenames = sorted([f for f in os.listdir(dr) if f.endswith('.cif')])

    mat = pd.DataFrame(index=filenames, columns=filenames, dtype=float)
    for f in filenames:
        mat.loc[f, f] = 0.0

    n = len(filenames)
    for i in range(n):
        f1 = filenames[i]
        id1 = f1[:-4]
        for j in range(i + 1, n):
            f2 = filenames[j]
            id2 = f2[:-4]

            if id1 not in res_order_map or id2 not in res_order_map:
                raise KeyError(f"Missing residue list for {id1} or {id2} in res_order_map")

            r = align_by_csv_order(
                f1, f2,
                res_order_map[id1],
                res_order_map[id2],
                dr,
                mode=mode,
                min_common_atoms=min_common_atoms
            )
            if decimals is not None:
                r = round(r, int(decimals))
            mat.loc[f1, f2] = r
            mat.loc[f2, f1] = r

    out = mat.reset_index().rename(columns={"index": "#"})
    return out


# functions to calculate and visualize clustering validation scores 
def compute_validation_scores(D1, Z, cutoff_range):
    rmsd_matrix = D1.values.astype(float)

    # NOTE: MDS can have randomness; set random_state for determinism
    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=0)
    coords = mds.fit_transform(rmsd_matrix)

    silhouette_scores = []
    calinski_scores = []
    db_scores = []
    valid_cutoffs = []

    for cutoff in cutoff_range:
        labels = fcluster(Z, t=cutoff, criterion='distance')
        n_clusters = len(set(labels))
        if n_clusters < 2 or n_clusters >= len(labels):
            continue

        valid_cutoffs.append(cutoff)
        silhouette_scores.append(silhouette_score(coords, labels))
        calinski_scores.append(calinski_harabasz_score(coords, labels))
        db_scores.append(davies_bouldin_score(coords, labels))

    validation_df = pd.DataFrame({
        'RMSD_cutoff': valid_cutoffs,
        'silhouette_score': silhouette_scores,
        'calinski_score': calinski_scores,
        'db_score': db_scores
    })
    return valid_cutoffs, silhouette_scores, calinski_scores, db_scores, validation_df

#this function will normalize the validation scores
def min_max_normalize(arr):
    arr = np.array(arr, dtype=np.float64)
    min_val = np.nanmin(arr)
    max_val = np.nanmax(arr)
    return (arr - min_val) / (max_val - min_val) if max_val > min_val else np.zeros_like(arr)
#this function will plot different validation scores against the RMSD cut-off ranges
def plot_validation_scores(cutoffs, sil_scores, ch_scores, db_scores, fstring):
    plt.figure(figsize=(12, 6))
    plt.rcParams['font.family'] = 'Arial'

    # Normalize each score to [0, 1] range
    sil_scores_norm = min_max_normalize(sil_scores)
    ch_scores_norm  = min_max_normalize(ch_scores)
    db_scores_norm  = min_max_normalize(db_scores)

    # Plotting normalized scores
    plt.plot(cutoffs, sil_scores_norm, marker='o', linestyle='-', label='Silhouette (norm)', color='royalblue')
    plt.plot(cutoffs, ch_scores_norm, marker='^', linestyle='-', label='Calinski-Harabasz (norm)', color='darkgreen')
    plt.plot(cutoffs, db_scores_norm, marker='s', linestyle='-', label='Davies-Bouldin (norm)', color='darkorange')

    # best cutoffs defined by different validation metrics
    best_cutoff_sil = cutoffs[int(np.nanargmax(sil_scores))]
    best_cutoff_ch  = cutoffs[int(np.nanargmax(ch_scores))]
    best_cutoff_db  = cutoffs[int(np.nanargmin(db_scores))]

    plt.axvline(x=best_cutoff_sil, color='royalblue', linestyle='-', linewidth=3, alpha=0.6,
                label=f'Best Silhouette = {best_cutoff_sil:.2f} Å')
    plt.axvline(x=best_cutoff_ch, color='darkgreen', linestyle='--', linewidth=2, alpha=0.6,
                label=f'Best Calinski = {best_cutoff_ch:.2f} Å')
    plt.axvline(x=best_cutoff_db, color='darkorange', linestyle='--', linewidth=1, alpha=0.6,
                label=f'Best Davies-Bouldin = {best_cutoff_db:.2f} Å')

    #adding labels for axis
    plt.xlabel("RMSD Cut-off (Å)", fontsize=18)
    plt.ylabel("Normalized Validation Score", fontsize=18)

    #modifying the legend box
    plt.legend(frameon=True)
    #turning on/off the grids in the plot
    plt.grid(False)

    ax = plt.gca()

    #modifying plot boarders
    ax.spines.top.set_linewidth(1)
    ax.spines.bottom.set_linewidth(1)
    ax.spines.left.set_linewidth(1)
    ax.spines.right.set_linewidth(1)

    #specifying decimal places for the axis tick labels
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    #modifying axis tick length and thickness
    ax.tick_params(axis='both', length=10, width=1)

    #modifying axis tick labels
    plt.xticks(rotation=0, fontsize=15)
    plt.yticks(rotation=0, fontsize=15)

    plt.tight_layout()
    plt.savefig(snakemake.output['val_score_plt'], format="png",
                bbox_inches="tight", dpi=2000)


# functions to calculate teh representative structure for each cluster

def load_structures_from_directory(directory, f_list):
    parser = MMCIFParser(QUIET=True)
    structures = []
    file_paths = [os.path.join(directory, f) for f in sorted(os.listdir(directory)) if f in f_list]
    for file_path in file_paths:
        structure = parser.get_structure(file_path, file_path)
        structures.append(structure)
    return structures

def get_strid_from_structure_obj(structure):
    sid = str(structure.id)
    sid = os.path.basename(sid)
    if sid.endswith(".cif"):
        sid = sid[:-4]
    return sid

def matched_atoms_for_same_strid(structure, csv_list, mode="backbone"):
    amap = atom_map_from_csv_order(structure, csv_list, mode=mode)
    keys = sorted(amap.keys())
    atoms = [amap[k] for k in keys]
    return keys, atoms

def align_structures_csv_order(structures, res_order_map,
                               mode="backbone", min_common_atoms=6):
    
    #Align all structures to the first structure (reference), using CSV-defined #correspondence.
    
    if not structures:
        raise ValueError("No structures provided.")

    ref = structures[0]
    ref_id = get_strid_from_structure_obj(ref)
    if ref_id not in res_order_map:
        raise KeyError(f"Missing residue order for reference {ref_id}")

    aligned = [ref]
    rmsds_to_ref = {ref_id: 0.0}

    for s in structures[1:]:
        sid = get_strid_from_structure_obj(s)
        if sid not in res_order_map:
            raise KeyError(f"Missing residue order for {sid}")

        atoms_ref, atoms_s = matched_atom_lists_from_csv(
            ref, s,
            res_order_map[ref_id],
            res_order_map[sid],
            mode=mode,
            min_common_atoms=min_common_atoms
        )

        sup = Superimposer()
        sup.set_atoms(atoms_ref, atoms_s)
        sup.apply(list(s.get_atoms()))  # move entire structure

        rmsds_to_ref[sid] = float(sup.rms)
        aligned.append(s)

    return aligned, rmsds_to_ref

def compute_average_structure_csv_order(ref_structure, aligned_structures,
                                       res_order_map, mode="backbone"):
    
    #Average coordinates for selected atoms (CSV-defined keys), returns avg_structure + #keys list.
    
    ref_id = get_strid_from_structure_obj(ref_structure)
    if ref_id not in res_order_map:
        raise KeyError(f"Missing residue order for reference {ref_id}")

    keys, ref_atoms = matched_atoms_for_same_strid(ref_structure, res_order_map[ref_id], mode=mode)
    if len(ref_atoms) == 0:
        raise ValueError("No atoms selected for averaging. Check mode/atom selection.")

    coords_sum = np.zeros((len(keys), 3), dtype=float)

    for s in aligned_structures:
        sid = get_strid_from_structure_obj(s)
        if sid not in res_order_map:
            raise KeyError(f"Missing residue order for {sid}")

        amap = atom_map_from_csv_order(s, res_order_map[sid], mode=mode)
        missing = [k for k in keys if k not in amap]
        if missing:
            raise ValueError(f"{sid} is missing {len(missing)} atoms needed for averaging.")

        coords = np.array([amap[k].get_coord() for k in keys], dtype=float)
        coords_sum += coords

    avg_coords = coords_sum / float(len(aligned_structures))

    avg_structure = ref_structure.copy()
    avg_map = atom_map_from_csv_order(avg_structure, res_order_map[ref_id], mode=mode)
    for i, k in enumerate(keys):
        avg_map[k].set_coord(avg_coords[i])

    return avg_structure, keys

def rmsd_to_average_csv_order_with_refid(avg_structure, ref_id, keys,
                                        structures, res_order_map,
                                        mode="backbone"):
    
    #RMSD(avg, structure_i) for each structure, using the same ordered keys.
    
    if ref_id not in res_order_map:
        raise KeyError(f"Missing residue order for reference {ref_id}")

    avg_map = atom_map_from_csv_order(avg_structure, res_order_map[ref_id], mode=mode)

    rmsd_dict = {}
    for s in structures:
        sid = get_strid_from_structure_obj(s)
        if sid not in res_order_map:
            raise KeyError(f"Missing residue order for {sid}")

        s_map = atom_map_from_csv_order(s, res_order_map[sid], mode=mode)
        missing = [k for k in keys if (k not in avg_map or k not in s_map)]
        if missing:
            raise ValueError(f"{sid} missing {len(missing)} atoms needed for RMSD-to-average.")

        avg_atoms = [avg_map[k] for k in keys]
        s_atoms   = [s_map[k] for k in keys]

        sup = Superimposer()
        sup.set_atoms(avg_atoms, s_atoms)
        rmsd_dict[sid] = float(sup.rms)

    return rmsd_dict


#################################### Analysis steps ####################################

def main():

    csvfile = snakemake.input["frag_info"]
    FSTRING = snakemake.input["frag_info"].split("/")[-1][:-14]
    val_mat = snakemake.config['val_score']
    leaf_lab = snakemake.config['leaf_labels']

    # ---------------- Step 1: import data ----------------
    DF1 = pd.read_csv(csvfile)
    DF1.replace(" ", "*", inplace=True)

    # fixing '5E54' becoming '5.00E+54' issue
    DF1 = fix_PDB_ID(DF1)
    
    #The following two lines will ensure that the datafrmae is sorted based on pdb, chain_A, resi_A_1 columns before generating unique structure IDs
    # This will also confirm that the order of the rows remains the same if a different order of the rows is provided in the input csv file
    sort_cols = ['pdb', 'chain_A', 'resi_A_1']
    DF1 = DF1.sort_values(sort_cols, kind='mergesort').reset_index(drop=True)

    # Make str_ID UNIQUE per row to prevent overwriting clipped CIFs
    DF1['str_ID'] = [
        f"{FSTRING}_{i:05d}_{p}_{c}_{r}"
        for i, (p, c, r) in enumerate(zip(DF1['pdb'], DF1['chain_A'], DF1['resi_A_1']))
    ]

    # ---------------- Step 2: generate clipped structures ----------------

    # Check whether the folder exists to write the clipped structures.
    # If it does not exist, create the folder.
    clip_path = snakemake.config['clipped_structures_dir'][:-1] + "_" + FSTRING + "/"
    if not os.path.isdir(clip_path):
        try:
            os.mkdir(clip_path)
        except FileExistsError:
            time.sleep(5.0)
            if not os.path.isdir(clip_path):
                os.mkdir(clip_path)

    # group columns like resi_A_1, resi_A_2
    suffix_groups = {}
    for col in DF1.columns:
        if len(col.split('_')) == 3:
            suffix = '_'.join(col.split("_")[1:])
            suffix_groups.setdefault(suffix, []).append(col)

    # Build and save residue-order map (CSV canonical order)
    res_order_map = {}

    for ind, pid in enumerate(DF1['pdb']):
        list_residue = []
        for suffix in suffix_groups:
            chain_col = 'chain_' + suffix.split('_')[0]
            res = DF1[chain_col][ind] + '_' + '_'.join([str(DF1[i][ind]) for i in suffix_groups[suffix]])
            list_residue.append(res)

        sid = DF1.loc[ind, 'str_ID']
        res_order_map[sid] = list_residue

        out_cif = sid + ".cif"
        print("Clipping:", out_cif)
        clip(pid, list_residue, 'N', out_cif, clip_path)

    with open(snakemake.output['res_order_map'], "w") as f:
        json.dump(res_order_map, f, indent=2)

    # ---------------- Step 3: all-vs-all RMSD ----------------

    # recommended: backbone-only correspondence (robust across sequence differences)
    RMSD_MODE = "backbone"
    # approx min atoms: 2 per residue
    any_sid = next(iter(res_order_map))
    min_common_atoms = 2 * len(res_order_map[any_sid])

    rmsd_df = align_all_against_all(
        clip_path,
        res_order_map=res_order_map,
        mode=RMSD_MODE,
        min_common_atoms=min_common_atoms,
        decimals=None
    )
    rmsd_df.to_csv(snakemake.output['rmsd'], index=False)

    # ---------------- Step 4: clustering + validation ----------------
    cl_names = [j.rstrip('.cif').replace('-', '_') for j in list(rmsd_df.columns)[1:]]

    D1 = rmsd_df.drop(['#'], axis=1)
    D1 = D1.astype(float)
    D1.index = cl_names
    D1.columns = cl_names

    # Condensed distance vector for linkage
    rmsd_condensed = squareform(D1.values)

    Z = linkage(rmsd_condensed, method='average')
    dend = dendrogram(Z, labels=D1.index, leaf_rotation=0, orientation='left', leaf_font_size=5, no_plot=True)

    min_rmsd1 = np.min(Z[:, 2])
    max_rmsd = np.max(Z[:, 2])
    pre_max = floor(max_rmsd)

    cutoffs_probe = np.linspace(min_rmsd1, pre_max, 50)
    min_rmsd = None
    for cutoff in cutoffs_probe:
        labels = fcluster(Z, t=cutoff, criterion='distance')
        cluster_sizes = pd.Series(labels).value_counts()
        if (cluster_sizes >= 4).any():
            min_rmsd = cutoff
            break

    print('#############################################################################################')
    print('Minimum RMSD cutoff with at least one cluster ≥ 4 members: ' + str(min_rmsd))
    print('#############################################################################################')

    cutoff_range = np.linspace(min_rmsd, pre_max, 50)

    cutoffs, sil_scores, ch_scores, db_scores, validation_df = compute_validation_scores(D1, Z, cutoff_range)
    plot_validation_scores(cutoffs, sil_scores, ch_scores, db_scores, FSTRING)
    validation_df.to_csv(snakemake.output['val_score_csv'], index=False)

    best_cutoff_sil = round(cutoffs[int(np.nanargmax(sil_scores))], 2)
    best_cutoff_ch = round(cutoffs[int(np.nanargmax(ch_scores))], 2)
    best_cutoff_db = round(cutoffs[int(np.nanargmin(db_scores))], 2)

    # ---------------- Dendrogram figure ----------------
    plt.figure(figsize=(3.5, 9))

    if val_mat == 'sil':
        dend = dendrogram(Z, color_threshold=best_cutoff_sil, labels=D1.index,
                          above_threshold_color='lightgray', leaf_rotation=0,
                          orientation='left', leaf_font_size=5, no_plot=True)
    elif val_mat == 'ch':
        dend = dendrogram(Z, color_threshold=best_cutoff_ch, labels=D1.index,
                          above_threshold_color='lightgray', leaf_rotation=0,
                          orientation='left', leaf_font_size=5, no_plot=True)
    elif val_mat == 'db':
        dend = dendrogram(Z, color_threshold=best_cutoff_db, labels=D1.index,
                          above_threshold_color='lightgray', leaf_rotation=0,
                          orientation='left', leaf_font_size=5, no_plot=True)
    else:
        raise ValueError("val_mat must be one of: sil, ch, or db")

    color_list_counter = dict(Counter(dend['color_list']))
    lcolor_list_counter = dict(Counter(dend['leaves_color_list']))

    small_clusters = []
    for c in lcolor_list_counter:
        if c != 'lightgray' and lcolor_list_counter[c] < 4:
            small_clusters.append(c)

    color_list_n = []
    leaves_color_list_n = []
    for cls in dend['color_list']:
        if cls in small_clusters or cls == 'lightgray':
            color_list_n.append('lightgray')
        else:
            color_list_n.append(cls)

    for ls in dend['leaves_color_list']:
        if ls in small_clusters or ls == 'lightgray':
            leaves_color_list_n.append('lightgray')
        else:
            leaves_color_list_n.append(ls)

    dend["color_list"] = color_list_n
    dend["leaves_color_list"] = leaves_color_list_n

    if not leaf_lab:
        dend1 = dendrogram(Z, color_threshold=0, labels=D1.index,
                           above_threshold_color='white', leaf_rotation=0,
                           orientation='left', leaf_font_size=5, no_labels=True)
    else:
        dend1 = dendrogram(Z, color_threshold=0, labels=D1.index,
                           above_threshold_color='white', leaf_rotation=0,
                           orientation='left', leaf_font_size=5)

    for i, dcoord in enumerate(dend['dcoord']):
        plt.plot(dcoord, dend['icoord'][i], color=dend['color_list'][i])

    for leaf, leaf_color in zip(plt.gca().get_yticklabels(), dend["leaves_color_list"]):
        leaf.set_color(leaf_color)

    for line in plt.gca().get_lines():
        line.set_linewidth(0.5)

    if val_mat == 'sil':
        plt.axvline(x=best_cutoff_sil, color='royalblue', linestyle='-', linewidth=3, alpha=0.6)
    elif val_mat == 'ch':
        plt.axvline(x=best_cutoff_ch, color='darkgreen', linestyle='--', linewidth=2, alpha=0.6)
    elif val_mat == 'db':
        plt.axvline(x=best_cutoff_db, color='darkorange', linestyle='--', linewidth=1, alpha=0.6)

    plt.xlabel('RMSD (Å)')
    plt.savefig(snakemake.output['dendrogram'], format="png", bbox_inches="tight", dpi=2000)

    # ---------------- Build cluster dictionary from dendrogram colors ----------------
    cl_dict = {}
    for idx, col in enumerate(dend['leaves_color_list']):
        if col == 'lightgray':
            continue
        cl_dict.setdefault(col, []).append(dend['ivl'][idx])

    # convert colors -> C1, C2, ...
    cl_dict = {f"{k[:1]}{i+1}": v for i, (k, v) in enumerate(cl_dict.items())}
    total = len(cl_dict)
    cl_dict = {f"C{total - i}": v for i, v in enumerate(cl_dict.values())}

    def get_key(value):
        for key, values in cl_dict.items():
            if value in values:
                return key
        return None

    DF1['cluster'] = DF1['str_ID'].apply(get_key)
    DF1['cluster'] = DF1['cluster'].fillna('unclustered')

    # ---------------- Step 5: Representative structures per cluster (UPDATED) ----------------
    DF1['rmsd_avg'] = np.nan
    DF1['rep_cls'] = np.nan

    for cls in cl_dict:
        str_ids = cl_dict[cls]  # these are str_IDs (no .cif)
        file_list = [sid + ".cif" for sid in str_ids]

        structures = load_structures_from_directory(clip_path, file_list)
        if len(structures) == 0:
            continue

        aligned_structures, _ = align_structures_csv_order(
            structures,
            res_order_map=res_order_map,
            mode=RMSD_MODE,
            min_common_atoms=min_common_atoms
        )

        ref_structure = aligned_structures[0]
        ref_id = get_strid_from_structure_obj(ref_structure)

        avg_str, keys = compute_average_structure_csv_order(
            ref_structure,
            aligned_structures,
            res_order_map=res_order_map,
            mode=RMSD_MODE
        )

        rmsds = rmsd_to_average_csv_order_with_refid(
            avg_structure=avg_str,
            ref_id=ref_id,
            keys=keys,
            structures=aligned_structures,
            res_order_map=res_order_map,
            mode=RMSD_MODE
        )

        # write to DF1
        DF1['rmsd_avg'] = DF1['rmsd_avg'].fillna(DF1['str_ID'].map(rmsds))

        rep_id = min(rmsds, key=rmsds.get)
        rmsd_rep = {sid: (1 if sid == rep_id else 0) for sid in rmsds}
        DF1['rep_cls'] = DF1['rep_cls'].fillna(DF1['str_ID'].map(rmsd_rep))

    DF1['rmsd_avg'] = DF1['rmsd_avg'].fillna('unclustered')
    DF1['rep_cls'] = DF1['rep_cls'].fillna('unclustered')
    DF1['rep_cls'] = DF1['rep_cls'].apply(lambda x: int(x) if isinstance(x, float) else x)

    # largest cluster
    if len(cl_dict) > 0:
        max_key = max(cl_dict, key=lambda k: len(cl_dict[k]))
    else:
        max_key = None

    DF1['rep_data'] = 0
    if max_key is not None:
        DF1['rep_data'] = ((DF1['cluster'] == max_key) & (DF1['rep_cls'] == 1)).astype(int)

    DF1.replace('*', '', inplace=True)
    DF1.to_csv(snakemake.output['align_clust_results'], index=False)

if __name__ == "__main__":
    main()
