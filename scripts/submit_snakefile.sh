#!/bin/bash
#SBATCH --job-name=snakemake_controller
#SBATCH --time=48:00:00
#SBATCH --mem=4G
#SBATCH --output=logs/snakemake_controller_%j.out
#SBATCH --partition=low

#run this scripts in running directory (usually one directory above /scripts)

module load conda/latest
cd /path/to/project/dir
source $(conda info --base)/etc/profile.d/conda.sh
conda activate snakemake_env

# debugging the logfile for conda init issues -- should version info in logfile(s) 
which snakemake
snakemake --version

echo "Conda executable: $CONDA_EXE"
echo "Conda level: $CONDA_SHLVL"

# use-envmodules is needed since I have some programs as modules, not in a conda environment
# we should move away from this approach soon, and load everything into conda environments
set -x
snakemake -s scripts/Snakefile.py --use-conda --rerun-incomplete --executor slurm --jobs 3 --printshellcmds --directory /path/to/project/dir --default-resources mem_mb=4096 runtime=600 slurm_partition=low
