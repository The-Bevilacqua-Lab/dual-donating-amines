# Dual-Donating Amines

This Snakemake workflow mines RNA-containing structures retrieved from the Protein Data Bank (PDB), identifies 
nucleobases that contain non-, single-, and dual-donating amines, and then performs a variety of analyses. To reduce 
the redundancy of structures, the workflow was designed to source information from csv files obtained from BGSU's 
[Representative Sets of RNA 3D Structures](https://rna.bgsu.edu/rna3dhub/nrlist). The open access paper 
corresponding to this workflow can be found at https://doi.org/10.1261/rna.080624.125.

## Running this Workflow

1. Install [Miniforge](https://github.com/conda-forge/miniforge) if you do not already have it installed.
2. Create a conda environment (named 'snakemake' in this example) with Snakemake installed by entering `mamba create -n 
snakemake bioconda::snakemake=9.9.0` into your terminal.
3. Activate the snakemake conda environment by entering `conda activate snakemake`.
4. If you intend to run this workflow on a cluster (recommended), enter `pip install 
snakemake-executor-plugin-cluster-generic`. This plugin will facilitate running Snakemake on the cluster.
5. Clone or download the dual-donating-amines repo and move to this folder within your terminal.
6. Consider whether you want to manually provide any mmCIF files for any specific PDB entry versions. If you do, 
   follow the instructions provided in the "PDB Entry Versions" section below.
7. If you are running the workflow on a cluster, follow the first substep below; otherwise, follow the 
   second substep below.
   1. I ran this workflow on Penn State's Roar Collab which uses the SLURM job manager. If the cluster you are using 
      also uses SLURM, review the contents of run.sh to ensure that the parameters are set to your liking, then 
      submit the script by entering `sbatch run.sh`. In particular, you may need to adjust the --jobs Snakemake 
      argument depending on how many cores are available. If you reduce this value, you may have to increase the 
      walltime specified on the fifth line. If you are dealing with a different job manager, you will need to modify 
      the run.sh script and how it is submitted to the queue.
   2. Enter `snakemake --sdm conda --cores 4` into your terminal. You can change the --cores argument depending on 
      how many cores you have available. I have not run the entire workflow in this way. Doing so with just a few 
      cores may take a long time to complete, possibly weeks.
8. After the run completes, the data will be available in the newly created results folder.

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
3. Move this file to resources/original_mmCIF_files/ within the dual-donating-amines folder.
4. Rename the file to 7k16.cif. On case-sensitive filesystems, it is important to use lowercase letters in the file 
   name.

The specific versions of the PDB entries that can be used to reproduce the results in the paper 
are provided in Supplemental_Table_S2.csv of the paper's 
[Supplemental Material](https://rnajournal.cshlp.org/content/32/1/21/suppl/DC1). In cases where the major version 
numbers currently available on [RCSB](https://www.rcsb.org/) are higher than what is 
reported in this CSV file, you may want to manually provide the mmCIF files with the lower major version numbers. Running this 
workflow without paying special attention to PDB entry versions should still produce results that are very close to 
what was published in the paper. However, providing mmCIF files with major version numbers that match what is listed 
in Supplemental_Table_S2.csv will improve alignment of the results. 
