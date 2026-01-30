#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --mem=5GB 
#SBATCH --time=48:00:00 
#SBATCH --account=open
#SBATCH --output=results/slurm/slurm-%j.out
#SBATCH --error=results/slurm/slurm-%j.out

snakemake --rerun-triggers mtime --sdm conda --jobs 500 --latency-wait 60 --executor cluster-generic --cluster-generic-submit-cmd 'sbatch --nodes=1 --ntasks=1 --mem=5GB --time=8:00:00 --account=open --output=results/slurm/slurm-%j.out --error=results/slurm/slurm-%j.out'
