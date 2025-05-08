

#################################### Importing libraries ####################################
import os
import shutil
import pandas as pd
import requests
import numpy as np
from Bio.PDB import MMCIFParser, Superimposer, PDBIO
from optparse import OptionParser
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from collections import Counter
from Bio.PDB import MMCIFParser, Superimposer, PDBIO
from optparse import OptionParser
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
from math import floor
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from sklearn.manifold import MDS
from optparse import OptionParser
from matplotlib.ticker import FormatStrFormatter

#################################### Functions ####################################

#fixing '5E54' becoming '5.00E+54' issue
def fix_PDB_ID(df):
    for i, j in enumerate(df['PDB_ID']):
        print (j)
        if len(str(j))>4:
            if '-' in j:
                X= ''.join(j.split('-')).upper()
                df.loc[i, 'PDB_ID']= X
            else:
                if '.' in j:
                    x= str(j)
                    x1= x.split('.')[0]
                    x2= x.split('+')[1]
                    X= x1+'E'+x2
                    print (X)
                    df.loc[i, 'PDB_ID']= X
                else:
                    x= str(j)
                    x1= x.split('+')[0]
                    x2= x.split('+')[1]
                    X= x1+ x2
                    print (X)
                    df.loc[i, 'PDB_ID']= X
                    
    return (df)

# functions to generate clipped structures

#function to extract structure file in cif format and convert as a pandas dataframe
def extract_cif(p):
    filename= p+'.cif'
    url= 'https://files.rcsb.org/view/%s' %filename 
    u = requests.get(url).content
    u1= str(u).split('#')
    if len(u1)==1:
        print ('STRUCTURE IS NOT AVAILABLE!')
        return '0'
    
    else:
        c_intrst= ['_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol',
           '_atom_site.label_atom_id', '_atom_site.label_alt_id',
           '_atom_site.label_comp_id', '_atom_site.label_asym_id',
           '_atom_site.label_entity_id', '_atom_site.label_seq_id',
           '_atom_site.pdbx_PDB_ins_code', '_atom_site.Cartn_x',
           '_atom_site.Cartn_y', '_atom_site.Cartn_z', '_atom_site.occupancy',
           '_atom_site.B_iso_or_equiv', '_atom_site.pdbx_formal_charge',
           '_atom_site.auth_seq_id', '_atom_site.auth_comp_id',
           '_atom_site.auth_asym_id', '_atom_site.auth_atom_id',
           '_atom_site.pdbx_PDB_model_num', '###']
        u2=[]
        for i, j in enumerate(u1):
            if c_intrst[0] in j:
                u2.append(j)

        if len(u2)==2:
            if len(u2[0])> len(u2[1]):
                u3= u2[0].replace('      ', '$')
            else:
                u3= u2[1].replace('      ', '$')
        else:
            u3= u2[0].replace('      ', '$')
        u4= u3.replace('     ', '$')
        u5= u4.replace('    ', '$')
        u6= u5.replace('   ', '$')
        u7= u6.replace('  ', '$')
        u8= u7.replace(' ', '$')
        u9= u8.replace('$$', '$')
        u10= u9.replace('\\\\', '').split('\\n')
        df = pd.DataFrame({'all_col':u10[23:]})

        df[c_intrst] = df['all_col'].str.split('$', expand=True)
        df= df.drop(['all_col', '###'], axis=1)

        return df


def clip(p, r, h, F):
    #p= pdb_ID, entries can be case-insensitive
    #r= list of residues 
    #each residue information in the list above will be formatted like: chain_ID.residue_ID
    #h= 'N' or 'Y', if 'N', that means no hydrogen coordinates should be included considering hydrogen coordinates are provided 
    #if 'Y', that means no hydrogen coordinates should be removed if any hydrogen coordinates are provided
    #F= desired output filename (with the extension of '.cif' or '.pdb')
    P= p.upper()
    c_ids= list(set([i.split('.')[0] for i in r]))
    print (c_ids)

    p1= extract_cif(P)

    if len(p1)>1:
        p1['temp_col'] = p1[['_atom_site.auth_asym_id', '_atom_site.auth_seq_id', '_atom_site.pdbx_PDB_ins_code']].fillna("").agg('_'.join, axis=1).str.strip()
        p1['temp_col'] = p1['temp_col'].str.replace('?', '*', regex=True)
        # instead of filtering out by chain ID, segment ID, residue index, residue ID one after another
        # temp column will have chain ID, residue index, and insertion code all together
        # each row will have their unique items for this column
        print (p1['temp_col'])
        p2= p1[p1['temp_col'].isin(r)]
        p2.index= np.arange(0, len(p2))

        #filtering out or not hydrogen atoms
        if h== 'N':
            p3= p2[p2['_atom_site.type_symbol']!= 'H']
            p3.index= np.arange(0, len(p3))
        else:
            p3= p2.copy()
        
        p3 = p3.drop(columns=['temp_col'])

        #cleaning name for the sugar atoms
        for i, j in enumerate(p3['_atom_site.label_atom_id']):
            if "\\'" in j:
                #print (j)
                #print (j)
                x= (j.replace("\\'", "'"))
                y= (x.replace('"', ''))
                p3.loc[i, '_atom_site.label_atom_id'] = y
                p3.loc[i, '_atom_site.auth_atom_id'] = y
        
        if F.split('.')[1]== 'pdb':
            #converting to pdb 
            columns_needed = ["_atom_site.id", "_atom_site.label_atom_id", "_atom_site.label_comp_id",
                          "_atom_site.auth_asym_id", "_atom_site.auth_seq_id", "_atom_site.pdbx_PDB_ins_code",
                          "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z",
                          "_atom_site.occupancy", "_atom_site.B_iso_or_equiv", "_atom_site.type_symbol"
                          ]
        
            # Ensure all required columns exist
            for col in columns_needed:
                if col not in p3.columns:
                    p3[col] = ""  # Fill missing columns with empty values

            # Open the output PDB file for writing
            with open(F, "w") as pdb_file:
                atom_serial_number = 1
                common_chain_ID= 'A' #temporary solution
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
                    f"{float(row['_atom_site.occupancy']):>6.2f}{float(row['_atom_site.B_iso_or_equiv']):>6.2f}          {row['_atom_site.type_symbol']:>2}\n")
                    atom_serial_number += 1
                
        else:
            #converting to cif
            #Creating dataframe for the column name
            c_intrst_1= ['_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol',
                '_atom_site.label_atom_id', '_atom_site.label_alt_id',
                '_atom_site.label_comp_id', '_atom_site.label_asym_id',
                '_atom_site.label_entity_id', '_atom_site.label_seq_id',
                '_atom_site.pdbx_PDB_ins_code', '_atom_site.Cartn_x',
                '_atom_site.Cartn_y', '_atom_site.Cartn_z', '_atom_site.occupancy',
                '_atom_site.B_iso_or_equiv', '_atom_site.pdbx_formal_charge',
                '_atom_site.auth_seq_id', '_atom_site.auth_comp_id',
                '_atom_site.auth_asym_id', '_atom_site.auth_atom_id',
                '_atom_site.pdbx_PDB_model_num']
            c_intrst_2= ['data_', '#','loop_', '_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol',
                '_atom_site.label_atom_id', '_atom_site.label_alt_id',
                '_atom_site.label_comp_id', '_atom_site.label_asym_id',
                '_atom_site.label_entity_id', '_atom_site.label_seq_id',
                '_atom_site.pdbx_PDB_ins_code', '_atom_site.Cartn_x',
                '_atom_site.Cartn_y', '_atom_site.Cartn_z', '_atom_site.occupancy',
                '_atom_site.B_iso_or_equiv', '_atom_site.pdbx_formal_charge',
                '_atom_site.auth_seq_id', '_atom_site.auth_comp_id',
                '_atom_site.auth_asym_id', '_atom_site.auth_atom_id',
                '_atom_site.pdbx_PDB_model_num']
    
            c_intrst_2[0]= 'data_'+P
    
            p4= pd.DataFrame(columns= c_intrst_1)
            p4['_atom_site.group_PDB']= c_intrst_2 #this p3 dataframe above p2 will be joined in p4. p3 will help pymol to visualize the file
    
            p5 = pd.concat([p4, p3], axis=0)

            print (p5.shape)

            p5.to_csv(F, header= False, sep=' ', index= False)

#functions for the alignment and distance calculation 

#this function will keep only the atoms required for the alignment
def req_atoms_alignment(structure):
    print (structure)
    atom_dicts= {'A': ["p", "C4'", "N9", "C6", "C2"], 'G': ["p", "C4'", "N9", "C6", "C2"], 'OMG': ["p", "C4'", "N9", "C6", "C2"], '2MG': ["p", "C4'", "N9", "C6", "C2"], 'C': ["p", "C4'", "N1", "C2", "C4"],'U': ["p", "C4'", "N1", "C2", "C4"], 'OMU': ["p", "C4'", "N1", "C2", "C4"], '5MC': ["p", "C4'", "N1", "C2", "C4"],'4AC': ["p", "C4'", "N1", "C2", "C4"], 'PSU': ["p", "C4'", "N1", "C2", "C4"], 'M7A': ["p", "C4'", "N9", "C6", "C2"], 'A2M': ["p", "C4'", "N9", "C6", "C2"], 'G2M': ["p", "C4'", "N9", "C6", "C2"], 'C2M': ["p", "C4'", "N1", "C2", "C4"],'U2M': ["p", "C4'", "N1", "C2", "C4"]}
    #residues_e= {res1-1: 1, res1: 2, res1+1: 3, res2-1: 4, res2: 5, res2+1: 6}
    #count=0
    good_atoms=[]
    for model in structure:
        for chain in model:
            for residue in chain:
                #if residues_e[residue.get_id()[1]] in co_ind:
                atoms= list(residue.get_atoms())
                    #atoms_noH=[]
                for atom in atoms:
                    #print ('----------------------------->.......look here')
                    #print (atom.get_name())
                    if atom.get_name() in atom_dicts[residue.get_resname()]:
                        good_atoms.append(atom)
                else:
                    pass
    #print (len(good_atoms))
    return good_atoms

#this function will align two structures and calculate RMSD between them
def align(str1, str2):
    parser = MMCIFParser(QUIET=True)
    superimposer = Superimposer()
    structure1 = parser.get_structure("RNA1", str1)
    structure2 = parser.get_structure("RNA2", str2)

    atom_dicts= {'A': ["p", "C4'", "N9", "C6", "C2"], 'G': ["p", "C4'", "N9", "C6", "C2"], 'C': ["p", "C4'", "N1", "C2", "C4"], 'U': ["p", "C4'", "N1", "C2", "C4"]}
    #reg3_8EV3_6.G105-6.U179.cif

    #str1_res1= int(str1.split('_')[2].split('-')[0][1:].rstrip('.cif'))
    #str1_res2= int(str1.split('_')[3][1:].rstrip('.cif'))

    #str2_res1= int(str2.split('_')[2].split('-')[0][1:].rstrip('.cif'))
    #str2_res2= int(str2.split('_')[3][1:].rstrip('.cif'))

    #find residues which will be considered to calculate RMSD
    #str1_good_res= who_is_missing(structure1, str1_res1, str1_res2) #good residues in structure1
    #str2_good_res= who_is_missing(structure2, str2_res1, str2_res2) #good residues in structure2

    #common_goods= list(set(str1_good_res).intersection(str2_good_res)) #common

    structure1_atoms= req_atoms_alignment(structure1)
    structure2_atoms= req_atoms_alignment(structure2)

    superimposer.set_atoms(structure1_atoms, structure2_atoms)
    rmsd = superimposer.rms


    return rmsd

#this function will calculate all-against-all RMSD between structures in a directory
def align_all_against_all(dr):
    home= os.getcwd()
    os.chdir(dr)
    clmns=['#']
    filenames=[]
    RMSDs={}
    for filename in os.listdir('.'):
        if filename.endswith('.cif'):
            clmns.append(filename)
            filenames.append(filename)
    print (filenames)
        
    t1= pd.DataFrame(columns = clmns)
    t1['#']= filenames

    def unique_pairs(items):
        pairs = []
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                pairs.append(str(items[i])+'*'+str(items[j]))
        return pairs
    unq_pairs= unique_pairs(filenames)

    for i, j in enumerate(unq_pairs):
        #unq_pairs_RMSD[j]=[]
        k1= j.split('*')[0]
        k2= j.split('*')[1]
        print ('------------------------------------------------------------>look here for issues')
        print (k1)
        print (k2)

        t1.loc[t1['#']== k1, k2]= round(align(k1, k2), 2)
        t1.loc[t1['#']== k2, k1]= round(align(k2, k1), 2)
        t1.loc[t1['#']== k1, k1]= 0
        t1.loc[t1['#']== k2, k2]= 0
        if round(align(k1, k2), 2) != round(align(k2, k1), 2):
            print ('<><><><><><><><><><><><><><><><><><><><><><><>')
            print (j)
            print (align(k1, k2))
            print (align(k2, k1))

    
    os.chdir(home)
    return t1


# functions to calculate and visualize clustering validation scores 

import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from sklearn.manifold import MDS
import matplotlib.pyplot as plt

def compute_validation_scores(D1, Z, cutoff_range):
    rmsd_matrix = D1.values
    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=0)
    coords = mds.fit_transform(rmsd_matrix)

    silhouette_scores = []
    calinski_scores = []
    db_scores = []
    valid_cutoffs = []

    for cutoff in cutoff_range:
        labels = fcluster(Z, t=cutoff, criterion='distance')
        n_clusters = len(set(labels))
        if n_clusters < 2 or n_clusters >= len(labels):  # skip single or all-different clusters
            silhouette_scores.append(np.nan)
            calinski_scores.append(np.nan)
            db_scores.append(np.nan)
            continue

        valid_cutoffs.append(cutoff)
        silhouette_scores.append(silhouette_score(coords, labels))
        calinski_scores.append(calinski_harabasz_score(coords, labels))
        db_scores.append(davies_bouldin_score(coords, labels))
        # Create a DataFrame from the results
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
    ch_scores_norm = min_max_normalize(ch_scores)
    db_scores_norm = min_max_normalize(db_scores)

    # Plotting normalized scores
    plt.plot(cutoffs, sil_scores_norm, marker='o', linestyle='-', label='Silhouette Score (normalized)', color='royalblue')
    plt.plot(cutoffs, ch_scores_norm, marker='^', linestyle='-', label='Calinski-Harabasz Index (normalized)', color='darkgreen')
    plt.plot(cutoffs, db_scores_norm, marker='s', linestyle='-', label='Davies-Bouldin Index (normalized)', color='darkorange')

    # best cutoffs defined by different validation metrics
    best_cutoff_sil = cutoffs[np.nanargmax(sil_scores)]
    best_cutoff_ch = cutoffs[np.nanargmax(ch_scores)]
    best_cutoff_db = cutoffs[np.nanargmin(db_scores)]

    plt.axvline(x=best_cutoff_sil, color='royalblue', linestyle='-', linewidth=6, alpha=0.4,
                label=f'Best Silhouette Cutoff = {best_cutoff_sil:.2f} Å')
    plt.axvline(x=best_cutoff_ch, color='darkgreen', linestyle='--', linewidth=4, alpha=0.4,
                label=f'Best Calinski-Harabasz Cutoff = {best_cutoff_ch:.2f} Å')
    plt.axvline(x=best_cutoff_db, color='darkorange', linestyle='--', linewidth=2, alpha=0.4,
                label=f'Best Davies-Bouldin Cutoff = {best_cutoff_db:.2f} Å')


    #adding labels for axis
    plt.xlabel("RMSD Cut-off (Å)", fontsize=18)
    plt.ylabel("Normalized Validation Score", fontsize=18)
    
    #modifying the legend box
    plt.legend(frameon= True) #add 'fontsize= n' if you want to change the fontsize 
    
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
    ax.tick_params(axis='both', length=10, width= 1)
    
    #modifying axis tick labels
    plt.xticks(rotation=0, fontsize=15)
    plt.yticks(rotation=0,fontsize=15)
    
    plt.tight_layout()
    plt.savefig(fstring+"_val_score_vs_RMSD_cutoff_normalized.png", format="png", bbox_inches="tight", dpi=2000)



# functions to calculate teh representative structure for each cluster
 
def load_structures_from_directory(directory, f_list):
    parser = MMCIFParser(QUIET=True)
    structures = []
    #file_paths = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.cif')]
    file_paths = [os.path.join(directory, f) for f in os.listdir(directory) if f in(f_list)]
    print (file_paths)

    for file_path in file_paths:
        structure = parser.get_structure(file_path, file_path)
        structures.append(structure)
    
    return structures


def align_structures(structures):
    #this function takes a list of structures as input
    #the first structure will be considered as the reference structure 
    #for all structures selected atoms (by req_atoms_alignment) will be considered for alignment
    
    #print (structures)
    ref_structure = structures[0]
    
    print ('this is reference structrure:')
    print (ref_structure)
    #ref_atoms = [atom for atom in ref_structure.get_atoms()]
    ref_atoms = req_atoms_alignment(ref_structure)
    
    super_imposers = []
    
    for structure in structures[1:]:
        filename = os.path.basename(str(structure)).rstrip('.cif>')
        print(filename)
        #atoms = [atom for atom in structure.get_atoms()]
        atoms = req_atoms_alignment(structure)
        super_imposer = Superimposer()
        super_imposer.set_atoms(ref_atoms, atoms)
        #super_imposer.apply(structure.get_atoms().get_coords())
        
        super_imposer.apply(structure.get_atoms())

        # Save the aligned version of structures
        ##io = io=MMCIFIO()
        ##io.set_structure(structure) 
        ##io.save(filename+"_aligned.cif")
        print ('Distance from the average: '+ str(round(super_imposer.rms, 2)))
        #print (super_imposer.rms)
        super_imposers.append(structure)
    
    return super_imposers


def compute_average_structure(structures, super_imposers):
    if not structures:
        raise ValueError("No structures provided.")
    
    ref_structure = structures[0]
    ref_atoms = req_atoms_alignment(ref_structure)
    num_atoms = len(ref_atoms)
    #print ('look here')
    #print (ref_atoms)
    
    coords_sum = np.zeros((num_atoms, 3))
    
    for i, structure in enumerate(structures):
        print ('working at the structure no. '+ str(i))
        print (structure)
        if i > 0:
            #super_imposers[i-1].apply(structure.get_atoms().get_coords())
            #super_imposers[i-1].get_atoms().get_coords()
            atoms = [atom for atom in req_atoms_alignment(super_imposers[i-1])]
            coords = np.array([atom.get_coord() for atom in atoms])
            #print (len(coords))
            coords_sum += coords
        else:
            atoms = [atom for atom in req_atoms_alignment(structure)]
            coords = np.array([atom.get_coord() for atom in atoms])
            #print (len(coords))
            coords_sum += coords
    
    avg_coords = coords_sum / len(structures)
    #print ((avg_coords[23]))
    
    avg_structure = ref_structure.copy()

    for i, atom in enumerate(req_atoms_alignment(avg_structure)):
        #print (i)
        #print (atom)
        #print ('----------------before')
        #print (atom.get_coord())
        atom.set_coord(avg_coords[i])
        #print ('----------------after')
        #print (atom.get_coord())

    return avg_structure


def align_avg(avg, strs):
    avg_atoms= req_atoms_alignment(avg)
    rmsd_dict={}
    for structure in strs:
        filename = os.path.basename(str(structure)).rstrip('.cif>')
        print(filename)
        
        atoms = req_atoms_alignment(structure)
        super_imposer = Superimposer()
        super_imposer.set_atoms(avg_atoms, atoms)
        
        rmsd_dict[filename]= super_imposer.rms
        
    return rmsd_dict

#################################### Analysis steps ####################################

# step 1: import data
optparser = OptionParser()
(options, args) = optparser.parse_args()
csvfile = args[0]
FSTRING= args[1]
DF1= pd.read_csv(csvfile) #this will be referred as the 'MAIN DATAFRAME'
DF1['str_ID'] = FSTRING+ '_'+ DF1[['pdb', 'chain_A', 'resi_A_1']].astype(str).agg('_'.join, axis=1)


# step 2: generate structures with clipped motifs
usr_dir= args[2]

work_path= os.path.join(usr_dir, FSTRING+'_identifying_data_rep_structures')
os.makedirs(work_path, exist_ok=True)

clip_path= os.path.join(work_path, FSTRING+'_all_clipped_structures') #all clipped structures will be stored here
os.makedirs(clip_path, exist_ok=True)

result_path= os.path.join(work_path, FSTRING+'_result_rep_structures') #result of alignment and clustering will be stored here
os.makedirs(result_path, exist_ok=True)

os.chdir(clip_path)

suffix_groups = {}
DF1.replace(" ","*", inplace=True)
for col in DF1.columns:
    if len(col.split('_'))==3:
        suffix = '_'.join(col.split("_")[1:]) # Extract last substring after "_"
        #print (suffix)
        suffix_groups.setdefault(suffix, []).append(col)

print (suffix_groups)

for ind, pid in enumerate(DF1['pdb']):
    print (ind)
    print (pid)
    list_residue=[]
    for suffix in suffix_groups:
        chain_col= 'chain_'+suffix.split('_')[0]
        #print (chain_col)
        #for i in suffix_groups[suffix]:
        #    print (df_r1[i][ind])
        #print (suffix_groups[suffix])
        #res= df_r1[chain_col][ind]+'.'+'.'.join([str(df_r1[i][ind]) for i in suffix_groups[suffix]])
        res= DF1[chain_col][ind]+'_'+'_'.join([str(DF1[i][ind]) for i in suffix_groups[suffix]])
        #print (res)
        ##res= [i for i in suffix_groups[suffix]]
        #list_residue.append(res.rstrip('*'))
        list_residue.append(res)
    print (list_residue)
    ###c_ids= list(set([i.split('.')[0] for i in list_residue])) #list of chain_IDs
    ###r_1= [res_ind.split('.')[1] for res_ind in list_residue if res_ind.split('.')[0]== c_ids[0]] #list of residue indeces

    F= FSTRING+'_'+pid+'_'+list_residue[0].rstrip('_*') + '.cif' #filename for the output/clipped structure files
    print (F)

    clip(pid, list_residue, 'N', F)

# step 3: align and calculate pair-wise distance
os.chdir(result_path)

rmsd_df= align_all_against_all(clip_path)
rmsd_df.to_csv(FSTRING+'_all_against_all_RMSD.csv', index= False)

#step 4: clustering
cl_names= [j.rstrip('.cif').replace('-', '_') for i,j in enumerate(list(rmsd_df.columns)[1:])]

D1= rmsd_df.drop(['#'], axis=1)
D1.index= cl_names
D1.columns= cl_names

# Convert DataFrame to squareform distance matrix
rmsd_matrix = squareform(D1.values)

# Perform hierarchical clustering
Z = linkage(rmsd_matrix, method='average')
Z[:, [0, 1]] = Z[:, [1, 0]] 
#plt.figure(figsize=(15, 40), dpi = 400)

#dend = dendrogram(Z, color_threshold= 3.5, labels=D1.index, above_threshold_color='lightgray', leaf_rotation=0, orientation='left', leaf_font_size=5)
dend = dendrogram(Z, labels=D1.index, leaf_rotation=0, orientation='left', leaf_font_size=5, no_plot=True)

# calculating the range for the cut-offs within which different cluster validation matrics will be calculated
min_rmsd = np.min(Z[:, 2]) #min value for the range
max_rmsd = np.max(Z[:, 2])
pre_max = floor(max_rmsd) # max value for the range
# generating RMSD cut-off range
cutoff_range = np.linspace(min_rmsd, pre_max, 50)

# calculating validation matrics 
cutoffs, sil_scores, ch_scores, db_scores, validation_df = compute_validation_scores(D1, Z, cutoff_range)

# visualizing validation matrics at different RMSD cut-off
plot_validation_scores(cutoffs, sil_scores, ch_scores, db_scores, FSTRING)
validation_df.to_csv(FSTRING+'_val_score_vs_RMSD_cutoff.csv', index= False)

# extracting best RMSD cut-off defined by different validation matrics
best_cutoff_sil = round(cutoffs[np.nanargmax(sil_scores)], 2)
best_cutoff_ch = round(cutoffs[np.nanargmax(ch_scores)], 2)
best_cutoff_db = round(cutoffs[np.nanargmax(db_scores)], 2)


# generating dendrogram with the RMSD cut-off best by the Slihouette score (default) or other matric defined by the user
plt.figure()
val_mat= args[3]
if val_mat in ['0', '1']:
    dend = dendrogram(Z, color_threshold= best_cutoff_sil, labels=D1.index, above_threshold_color='lightgray', leaf_rotation=0, orientation='left', leaf_font_size=5, no_plot= True)
elif val_mat == '2':
    dend = dendrogram(Z, color_threshold= best_cutoff_ch, labels=D1.index, above_threshold_color='lightgray', leaf_rotation=0, orientation='left', leaf_font_size=5, no_plot= True)
elif val_mat == '3':
    dend = dendrogram(Z, color_threshold= best_cutoff_db, labels=D1.index, above_threshold_color='lightgray', leaf_rotation=0, orientation='left', leaf_font_size=5, no_plot= True)

color_list_counter=  dict(Counter(dend['color_list']))
lcolor_list_counter=  dict(Counter(dend['leaves_color_list']))

small_clusters=[] #clusters with less than 4 members will be stored here

for i in lcolor_list_counter:
    if i != 'lightgray' and lcolor_list_counter[i]<4:
        small_clusters.append(i)
print ('------------------------------------')
print (small_clusters)
        
color_list_n= [] #updated color list for the branches
leaves_color_list_n= [] #updated color list for the branch leaves

for i, cls in enumerate(dend['color_list']):
    if cls in small_clusters:
        color_list_n.append('lightgray')
    else:
        color_list_n.append(cls)

for i, ls in enumerate(dend['leaves_color_list']):
    if ls in small_clusters:
      leaves_color_list_n.append('lightgray')
    else:
        leaves_color_list_n.append(ls)

dend["color_list"]= color_list_n #updated branch color list is assigned to the original dendrogram
dend["leaves_color_list"]= leaves_color_list_n #updated leaves color list is assigned to the original dendrogram

leaf_lab= args[4]
if leaf_lab== '0': #this means user does not want leaves labeling 
    dend1 = dendrogram(Z, color_threshold= best_cutoff_sil, labels=D1.index, above_threshold_color='lightgray', leaf_rotation=0, orientation='left', leaf_font_size=5, no_labels=True)
elif leaf_lab== '1': #this means user wants the leaves labeling
    dend1 = dendrogram(Z, color_threshold= best_cutoff_sil, labels=D1.index, above_threshold_color='lightgray', leaf_rotation=0, orientation='left', leaf_font_size=5)

#assigning branch color
for i, dcoord in enumerate(dend['dcoord']):
    plt.plot(dcoord, dend['icoord'][i], color=dend['color_list'][i])

#assigning branch leaves color
for leaf, leaf_color in zip(plt.gca().get_yticklabels(), dend["leaves_color_list"]):
    leaf.set_color(leaf_color)

#modifying the branch thickness in the dendrogram
for line in plt.gca().get_lines():
    line.set_linewidth(1) 

if val_mat in ['0', '1']:
    plt.axvline(x=best_cutoff_sil, color='green', linestyle='--', alpha= 0.8,  linewidth=2)
elif val_mat=='2':
    plt.axvline(x=best_cutoff_ch, color='green', linestyle='--', alpha= 0.8,  linewidth=2)
elif val_mat=='3':
    plt.axvline(x=best_cutoff_db, color='green', linestyle='--', alpha= 0.8,  linewidth=2)

plt.savefig(FSTRING+"_all_structure_dendrogram.png", format="png", bbox_inches="tight", dpi = 2000)
#need to turn off the leaves labeling

cl_dict={}

for cls1, cls2 in enumerate(dend['leaves_color_list']):
    if cls2== 'lightgray':
        pass
    else:
        if cls2 in cl_dict:
            cl_dict[cls2].append(dend['ivl'][cls1])
        else:
            cl_dict[cls2]= []
            cl_dict[cls2].append(dend['ivl'][cls1]) 
cl_dict= {f"{k[:1]}{i+1}": v for i, (k, v) in enumerate(cl_dict.items()) if k != 'lightgray'}


# adding cluster information to the 'MAIN DATAFRAME'
def get_key(value):
    for key, values in cl_dict.items():
        if value in values:
            return key      
DF1['cluster'] = DF1['str_ID'].apply(get_key)
DF1['cluster'] = DF1['cluster'].fillna('unclustered')

# step 5: identifying representative structures for each cluster

DF1['rmsd_avg']= np.nan
DF1['rep_cls']= np.nan

for cls in cl_dict:
    str_list= [i+'.cif' for i in cl_dict[cls]]
    print (str_list)
    
    print ('STARTING STEP 1---------------------------------------------->>>>>>>>>>>>>')
    structures1= load_structures_from_directory(clip_path, str_list) #loading required structures

    print ('STARTING STEP 2---------------------------------------------->>>>>>>>>>>>>')
    structures2= align_structures(structures1) #aligning structures with the reference structure

    print ('STARTING STEP 3---------------------------------------------->>>>>>>>>>>>>')
    avg_str= compute_average_structure(structures1, structures2) #calculating average structure

    print ('STARTING STEP 4---------------------------------------------->>>>>>>>>>>>>')
    rmsds= align_avg(avg_str, structures1) #calculating distance from the average structure
    
    print (rmsds)
    
    #small_df['rmsd_avg'] = small_df['str_ID'].map(rmsds)
    DF1['rmsd_avg'] = DF1['rmsd_avg'].fillna(DF1['str_ID'].map(rmsds))
    #small_df['rmsd_avg'] = small_df['rmsd_avg'].fillna('unclustered')
    
    min_key = min(rmsds, key=rmsds.get)
    print ('-----------------------------------')
    print (min_key)
    
    rmsd_rep= {}
    
    for strs in rmsds:
        if strs== min_key:
            rmsd_rep[strs]= int(1)
        else:
            rmsd_rep[strs]= int(0)
    print (rmsd_rep)      
    #small_df['rep_cls'] = small_df['str_ID'].map(rmsd_rep)
    DF1['rep_cls'] = DF1['rep_cls'].fillna(DF1['str_ID'].map(rmsd_rep))
    #small_df['rep_cls'] = small_df['rep_cls'].fillna('unclustered')
    
    #small_df['rep_cls'] = small_df['rep_cls'].apply(lambda x: int(x) if isinstance(x, float) else x)
    
DF1['rmsd_avg'] = DF1['rmsd_avg'].fillna('unclustered')
DF1['rep_cls'] = DF1['rep_cls'].fillna('unclustered')
DF1['rep_cls'] = DF1['rep_cls'].apply(lambda x: int(x) if isinstance(x, float) else x)

#find the cluster with max number of structures
max_key = max(cl_dict, key=lambda k: len(cl_dict[k]))

# identifying data representative structure which is the representative structure for the major cluster
# major cluster is the cluster with the highest number of members
DF1['rep_data'] = ((DF1['cluster'] == max_key) & (DF1['rep_cls'] == 1)).astype(int)
DF1.to_csv(FSTRING+'_alignment_n_clustering_results.csv', index= False)

#  storing the data representative structure to the result_path
data_rep_str= DF1[DF1['rep_data']==1]['str_ID'].to_list()[0]

src = clip_path+ '/'+ data_rep_str+'.cif'
dst = result_path+ '/'+ data_rep_str+ '.cif'

shutil.copy(src, dst)