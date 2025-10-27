# Lemay Lab MAG reconstruction pipeline (Snakemake)

A reproducible Snakemake workflow for MAG reconstruction from FL100 cleaned reads and contigs. This repository was originally developed by Nithya K Kumar for metagenomic read preprocessing, and was modified by Margaret A. Hilliard using parts of the [SnakeMAGs](https://github.com/Nachida08/SnakeMAGs) framework. The original analysis repository can be found [here](https://github.com/margarethilliard/VB12-analysis). Metagenomic sequence pre-processing steps were run previously to generate input files, and that repository can be found [here](https://github.com/dglemay/ARG_metagenome).

---

## Overview
This pipeline streamlines the reconstruction of MAGs from cleaned reads and contigs. It performs binning, quality control and filtering, taxonomic classification and an abundance calculation in a reproducible and modular framework.

**Core features:**
- Automated workflow management using **Snakemake**
- Reproducible environments using **Conda**
- Modular structure for easy extension
- Example configuration for quick setup

---

##  Workflow summary

| Step | Tool | Description |
|------|------|--------------|
| 1. Mapping | BWA | Map reads to contigs |
| 2. Sort/compress | Samtools | Sort and compress SAM to BAM files |
| 3. Bin | Metabat2 | Generate read depth file and bin contigs |
| 4. QC | CheckM | Assess bin completion, contamination, and filters accordingly |
| 5. Taxonomic profiling | GTDB-tk | Taxonomic classification based on the Genome Database Taxonomy |
| 6. MAG abundances | CoverM | Generate read coverage and relative abundance information |

---

## Clone the repository
```bash
git clone https://github.com/margarethilliard/lemaylab_MAG_reconstruction
cd lemaylab_MAG_reconstruction
```


 Usage instructions
 ------------------
1.   Navigate to project root directory:
```bash
cd /path/to/project_root
```
 
2. Copy config/config.yaml and edit paths for your system. **This should be the only file you need to edit, you should not need to edit the Snakefile**

3. Create tab-separated sample sheet (sample_sheet.txt) with columns: sample_name, contig_path, r1_path, r2_path
```bash
sample_name	contig_path	r1_path	r2_path
sample01	/path/to/sample01_assembled/final.contigs.fa	/path/to/sample01_R1_dup.fastq.gz	/path/to/sample01_R2_dup.fastq.gz
sample02	/path/to/sample02_assembled/final.contigs.fa.gz	/path/to/sample02_R1_dup.fastq.gz	/path/to/sample02_R2_dup.fastq.gz
sample03	/path/to/sample03_assembled/final.contigs.fa.gz	/path/to/sample03_R1_dup.fastq.gz	/path/to/sample03_R2_dup.fastq.gz
```

4. Download the latest release of GTDB-Tk, which requires ~66G+ of external data (GTDB) that need to be downloaded and unarchived. 
```bash
# Download 
wget https://data.gtdb.ecogenomic.org/releases/release214/214.1/auxillary_files/gtdbtk_r214_data.tar.gz
# Decompress
tar -xzvf *tar.gz
# This will create a folder called release214
```

5. Load snakemake v9.11.4 into a conda environment (if necessary):
```bash
eval "$(mamba shell hook --shell bash)"
mamba create -n snakemake_env -c conda-forge -c bioconda snakemake=9.11.4
conda activate snakemake_env
```

 6. Install slurm executor plugin for snakemake v8+ (only needs to be done once):
 ```bash
pip install snakemake-executor-plugin-slurm
```

 7. Quick check to make sure there are no errors (dry run):
```bash
snakemake -s scripts/Snakefile --configfile config.yaml -n
```

8. Run the pipeline using one of these methods (meant for using HPC with SLURM scheduler):
* METHOD A - Submit via sbatch script (recommended):
    ```
    sbatch scripts/submit_snakefile.sh
    ````

* METOD B - Run in terminal directly. 
    ```bash
    snakemake -s scripts/Snakefile --configfile config.yaml --executor slurm --jobs 20 --use-conda \
            --default-resources slurm_account=GROUPNAME mem_mb=4096 runtime=600
    ```

9. Monitor progress:
```bash
tail -f logs/snakemake_<jobid>.out
```

Required directory structure:
-----------------------
```
# project_root/
# │
# ├── config/
# │   └── config.yaml               [REQUIRED - edit with your paths]
# │
# ├── scripts/
# │   ├── Snakefile                 [this file]
# │   └── submit_snakefile.sh       [submit this file with sbatch]
# │
# ├── envs/
# │   ├── BWA.yaml              
# │   ├── CHECKM.yaml   
# │   ├── COVERM.yaml                
# │   ├── GTDBTK.yaml 
# │   ├── METABAT2.yaml        
# │	  └── SAMTOOLS.yaml 
# │	
# ├── sample_sheet.txt             [REQUIRED - tab-separated file with sample info]
# │
# ├──fastq_files/                  [your input files from fastuniq]
# │   ├── sample01_R1.fastq.gz
# │   ├── sample01_R2.fastq.gz
# │   ├── sample02_R1.fastq.gz
# │   ├── sample02_R2.fastq.gz
# │   ├── sample03_R1.fastq.gz
# │   └── sample03_R2.fastq.gz	
# │
# └── contig_files/                [your input files from MEGAHIT]
#     ├── sample1_final.contigs.fa.gz
#     ├── sample2_final.contigs.fa.gz
#     └── sample3_final.contigs.fa.gz
```
