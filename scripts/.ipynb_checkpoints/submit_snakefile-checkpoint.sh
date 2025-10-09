#!/bin/bash
#SBATCH --job-name=snakemake_controller_igaseq
#SBATCH --time=48:00:00
#SBATCH --mem=4G
#SBATCH --output=logs/snakemake_controller_%j.out
##SBATCH --account=dglemaygrp
#SBATCH --partition=low

#run this scripts in running directory (usually one directory above /scripts)

source ~/.bashrc
conda activate snakemake_env

# use-envmodules is needed since I have some programs as modules, not in a conda environment
# we should move away from this approach soon, and load everything into conda environments 
snakemake -s scripts/Snakefile_step1_3_FINAL.py \
    --executor slurm \
    --jobs 20 \
    --use-conda \
    --rerun-incomplete \
    --printshellcmds \
    --default-resources mem_mb=4096 runtime=600 slurm_partition=low