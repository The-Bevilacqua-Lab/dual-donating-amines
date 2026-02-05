# Dual-Donating Amines

This Snakemake workflow mines RNA-containing structures retrieved from the Protein Data Bank (PDB), identifies 
nucleobases that contain non-, single-, and dual-donating amines, and then performs a variety of analyses. To reduce 
the redundancy of structures, the workflow was designed to source information from csv files obtained from BGSU's 
[Representative Sets of RNA 3D Structures](https://rna.bgsu.edu/rna3dhub/nrlist). The open access paper 
corresponding to this workflow can be found at https://doi.org/10.1261/rna.080624.125, and some additional related 
materials can be found on [figshare](https://figshare.com/s/1da9267f13b4be217871?file=61592395).

## Running this Workflow

Note: This workflow works best on case-sensitive filesystems. It should still run on case-insensitive filesystems 
(e.g., APFS), but the results may differ some due to file overwriting.

1. Install NaTorsion, available from [AMIGOSIII](https://github.com/pylelab/AMIGOSIII) (v1.2-alpha). 
2. Install [Miniforge](https://github.com/conda-forge/miniforge) if you do not already have it installed. Release 
   25.11.0-1 worked for me, although others will likely work as well.
3. Create a conda environment (named 'snakemake' in this example) with Snakemake installed by entering `mamba create -n 
   snakemake bioconda::snakemake=9.9.0` into your terminal.
4. Activate the snakemake conda environment by entering `conda activate snakemake`.
5. If you intend to run this workflow on a cluster (recommended), enter `pip install 
   snakemake-executor-plugin-cluster-generic`. This plugin will facilitate running Snakemake on the cluster.
6. Clone or download the dual-donating-amines repo and move to this folder within your terminal.
7. Consider whether you want to manually provide any mmCIF files for any specific PDB entry versions. If you do, 
   follow the instructions provided in the "PDB Entry Versions" section below.
8. If you are running the workflow on a cluster, follow the first substep below; otherwise, follow the 
   second substep below.
   1. I ran this workflow on Penn State's Roar Collab which uses the SLURM job manager. If the cluster you are using 
      also uses SLURM, review the contents of run.sh to ensure that the parameters are set to your liking, then 
      submit the script by entering `sbatch run.sh`. In particular, you may need to adjust the --jobs Snakemake 
      argument depending on how many cores are available. If you reduce this value, you may have to increase the 
      walltime specified on the fifth line. If you are dealing with a different job manager, you will need to modify 
      the run.sh script and how it is submitted to the queue.
   2. Enter `snakemake --rerun-triggers mtime --sdm conda --cores 4` into your terminal. You can change the --cores 
      argument depending on how many cores you have available. I have not run the entire workflow in this way. Doing 
      so with just a few cores may take a long time to finish.
9. After the run completes, the data will be available in the newly created results folder.

## PDB Entry Versions

By default, the workflow considers the latest versions of PDB entries; however, specific mmCIF files that correspond 
to different versions of PDB entries can also be provided. I will walk through an example with PDB ID 7K16. At the 
time of writing, the 
latest version available for this PDB entry is 2.0, dated May 21, 2025. Meanwhile, version 1.2, dated October 18, 
2023, is available through [RCSB](https://www.rcsb.org/). By performing the steps below, the workflow will consider 
the older version of this PDB entry.

1. On the [version webpage](https://www.rcsb.org/versions/7K16) for PDB ID 7K16, download the mmCIF file 
   corresponding to version 1.2.
2. Decompress the downloaded file in your terminal by entering `gunzip pdb_00007k16_xyz_v1-2.cif.gz`.
3. Create a folder named "original_mmCIF_files" within dual-donating-amines/resources/, then move the CIF file to 
   this new folder.
4. Rename the file to 7k16.cif. On case-sensitive filesystems, it is important to use lowercase letters in the file 
   name.

The specific versions of the PDB entries that can be used to reproduce the results in the paper are provided in 
Supplemental_Table_S2.csv of the paper's 
[Supplemental Material](https://rnajournal.cshlp.org/content/32/1/21/suppl/DC1). In cases where the major version 
numbers currently available on [RCSB](https://www.rcsb.org/) are higher than what is reported in this CSV file, you may 
want to manually provide the mmCIF files with the lower major version numbers. Running this workflow without paying 
special attention to PDB entry versions should still produce results that are very close to what was published in the 
paper. However, providing mmCIF files with major version numbers that match what is listed in Supplemental_Table_S2.csv 
will improve alignment of the results. 

## Comments on Clustering

### Original Method

The clustering method used to generate the results presented in the paper is sensitive to the ordering of rows in the 
fragment information CSV file that is used as the input to id_rep_struct.py. To improve reproducibility, code was added 
to create_plots.R, which applies a custom row ordering for the Location 1 and Location 2 data. Whether or not this 
code is run is controlled by the custom_order variable (defaults to true) in config.yaml.

Furthermore, during the analysis that generated the plots for the paper, the PDB ID for one of the structures (PDB ID 
9E71) was represented in scientific notation within the Location 2 data ("9e+71"). In other words, "9E71" was 
automatically stored as float64 in Python and as double in R. The workflow was since then revised such that PDB IDs 
are now stored as str in Python and as character in R.  

We later noticed that our original clustering method skipped PDB ID 9E71 because it was represented in scientific 
notation during the earlier analysis. To improve reproducibility, code was added to create_plots.R, pt_loc_info.py, 
and id_rep_struct.py to convert PDB ID 9E71 into scientific notation and to handle this structure similar to the 
previous analysis. Whether PDB ID 9E71 is converted into "9e+71" is controlled by the 9E71_sci_not variable 
(defaults to true) in config.yaml.

In sum, to generate clustering analysis plots that closely match what is depicted in the paper, leave the variables 
custom_order and 9E71_sci_not in config.yaml set to their default values of true.

### Improved Method

After realizing that the original clustering method is sensitive to row ordering, we opted to provide 
id_rep_struct_v2.py, which includes the following improvements:

- **Row Sorting:** The rows of the input data are automatically sorted such that differences in the original row 
  order will not lead to differences in the resulting plots.
- **PDB IDs and Scientific Notation:** A function named fix_PDB_ID is called to convert any PDB IDs that are 
  represented in scientific notation back to their typical four-character format. 
- **Atom-Atom Mapping:** Prior to performing the alignment for each pair of fragments, a mapping between the 
  coarse-grained atoms of the two fragments is performed. For instance, the C2 coarse-grained atom of a particular 
  adenine in the first fragment will be matched with the corresponding coarse-grained atom and residue from the second 
  fragment before aligning the two fragments to one another. This results in lower RMSD values.
- **No Rounding:** The original approach rounded the RMSD values to the hundredths place before clustering those 
  values. This introduced the possibility for artificial ties between different pairs of fragments. The RMSD values 
  are no longer rounded prior to clustering. 

As was the case with the original method, application of this improved method to the Location 1 and Location 2 data 
reveals representative structures that exhibit S-motifs. The plots that result from this new clustering appraoch are 
provided in [figshare](https://figshare.com/s/1da9267f13b4be217871?file=61592395). This improved clustering method 
can be used in this workflow by setting the original_clustering variable in config.yaml to false.
