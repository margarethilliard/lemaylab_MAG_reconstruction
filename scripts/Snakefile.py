# ------------------
# Usage instructions
# ------------------
# 1. Navigate to project root directory:
# 	 mkdir /quobyte/dglemaygrp/mhilliard/FL100/lemaylab_MAG_reconstruction && cd /quobyte/dglemaygrp/mhilliard/FL100/lemaylab_MAG_reconstruction
#    git clone https://github.com/margarethilliard/lemaylab_MAG_reconstruction.git .
# 
# 2. Edit config/config.yaml paths for your system
#
# 3. Copy fastq reads and contigs into fastq_files/ and contigs/ directories 
#
# 3. Create sample sheet (sample_sheet.txt) with columns: sample_name, contigs_path, r1_path, r2_path
#    See example at bottom of this file
#
# 4. Load snakemake v9.11.4 into a conda environment (if necessary):
#	 module load conda/latest
#	 eval "$(mamba shell hook --shell bash)"
#	 mamba create -n snakemake_env -c conda-forge -c bioconda snakemake=9.11.4
#	 mamba activate snakemake_env	
#
# 5. Install slurm executor plugin for snakemake v8+ (only needs to be done once):
#    pip install snakemake-executor-plugin-slurm
#
# 6. Run the pipeline using one of these methods:
#
#    METHOD A - Submit via sbatch script (recommended):
#    sbatch scripts/submit_snakefile.sh
#
#    METHOD B - Run directly with snakemake:
#    snakemake -s scripts/Snakefile --configfile config.yaml --executor slurm --jobs 20 --use-conda \
#        --default-resources slurm_account=GROUPNAME mem_mb=4096 runtime=600
#
# 7. Monitor progress:
#    tail -f logs/snakemake_<jobid>.out
#
########################################
# Required directory structure: 
# 
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
# │   ├── sample1_R1.fastq.gz
# │   ├── sample1_R2.fastq.gz
# │   ├── sample2_R1.fastq.gz
# │   ├── sample2_R2.fastq.gz
# │   ├── sample3_R1.fastq.gz
# │   └── sample3_R2.fastq.gz	
# │
# └── contig_files/                [your input files from MEGAHIT]
#     ├── sample1_final.contigs.fa.gz
#     ├── sample2_final.contigs.fa.gz
#     └── sample3_final.contigs.fa.gz
# 
# All output directories will be created automatically by Snakemake
#
########################################
# SAMPLE SHEET FORMAT (sample_sheet.txt):
# Tab-separated file with these columns:
#sample_name	contig_path	r1_path	r2_path
#5001	/path/to/5001_assembled/5001.final.contigs.fa	/path/to/5001_R1_dup.fastq.gz	/path/to/5001_R2_dup.fastq.gz
#5002	/path/to/5002_assembled/5002.final.contigs.fa.gz	/path/to/5002_R1_dup.fastq.gz	/path/to/5002_R2_dup.fastq.gz
#5006	/path/to/5006_assembled/5006.final.contigs.fa.gz	/path/to/5006_R1_dup.fastq.gz	/path/to/5006_R2_dup.fastq.gz
#
########################################
# Snakefile for FL100 MAG reconstruction
# Steps: BWA → Samtools → MetaBAT2 → CheckM
#                       └──→ Filtered bins (only QC-passed MAGs)
#                                ├──→ GTDB-Tk
#                                └──→ CoverM
########################################

import os
import pandas as pd

# ---- CONFIGURATION ----
configfile: "config/config.yaml"

PROJECT_ROOT = config["project_root"]

# Load sample sheet
samples_df = pd.read_csv(os.path.join(PROJECT_ROOT, config["sample_sheet"]), sep="\t")
# Force sample name to string 
samples_df["sample_name"] = samples_df["sample_name"].astype(str)
SAMPLES = samples_df["sample_name"].tolist()

SAMPLE_R1 = dict(zip(samples_df["sample_name"], [os.path.join(PROJECT_ROOT, p) for p in samples_df["r1_path"]]))
SAMPLE_R2 = dict(zip(samples_df["sample_name"], [os.path.join(PROJECT_ROOT, p) for p in samples_df["r2_path"]]))
SAMPLE_CONTIGS = dict(zip(samples_df["sample_name"], [os.path.join(PROJECT_ROOT, p) for p in samples_df["contigs_path"]]))

# Constants
GB = 1024

# ---- RULE ALL ----
rule all:
    input:
        expand(os.path.join(PROJECT_ROOT, "results/coverm/{sample}/abundance.tsv"), sample=SAMPLES),
        os.path.join(PROJECT_ROOT, "results/checkm/QC_summary.tsv")

########################################
# 1. BWA: Map reads to contigs
########################################
rule bwa_map:
    input:
        contigs = lambda wc: SAMPLE_CONTIGS[wc.sample],
        r1 = lambda wc: SAMPLE_R1[wc.sample],
        r2 = lambda wc: SAMPLE_R2[wc.sample]
    output:
        sam = os.path.join(PROJECT_ROOT, "results/bwa/{sample}/{sample}.sam")
    conda:
        lambda wc: os.path.join(PROJECT_ROOT, config["conda_envs"]["bwa"])
    threads: config["resources"]["bwa"]["threads"]
    resources:
        mem_mb = config["resources"]["bwa"]["mem_mb"],
        runtime = config["resources"]["bwa"]["runtime"]
    shell:
        f"""
        mkdir -p {PROJECT_ROOT}/results/bwa/{{wildcards.sample}}
        bwa index {{input.contigs}}
        bwa mem -t {{threads}} {{input.contigs}} {{input.r1}} {{input.r2}} > {{output.sam}}
        """

########################################
# 2. Samtools: Sort and compress SAM to BAM
########################################
rule samtools_sort:
    input:
        lambda wc: os.path.join(PROJECT_ROOT, f"results/bwa/{wc.sample}/{wc.sample}.sam")
    output:
        bam = os.path.join(PROJECT_ROOT, "results/samtools/{sample}/{sample}.sorted.bam")
    conda:
        lambda wc: os.path.join(PROJECT_ROOT, config["conda_envs"]["samtools"])
    threads: config["resources"]["samtools"]["threads"]
    resources:
        mem_mb = config["resources"]["samtools"]["mem_mb"],
        runtime = config["resources"]["samtools"]["runtime"]
    shell:
        f"""
        mkdir -p {PROJECT_ROOT}/results/samtools/{{wildcards.sample}}
        samtools sort -@ {{threads}} -o {{output.bam}} {{input}}
        samtools index {{output.bam}}
        """

########################################
# 3. MetaBAT2: Binning
########################################
rule metabat2:
    input:
        contigs = lambda wc: SAMPLE_CONTIGS[wc.sample],
        bam = lambda wc: os.path.join(PROJECT_ROOT, f"results/samtools/{wc.sample}/{wc.sample}.sorted.bam")
    output:
        bins_dir = directory(os.path.join(PROJECT_ROOT, "results/metabat2/{sample}/bins"))
    conda:
        lambda wc: os.path.join(PROJECT_ROOT, config["conda_envs"]["metabat2"])
    threads: config["resources"]["metabat2"]["threads"]
    resources:
        mem_mb = config["resources"]["metabat2"]["mem_mb"],
        runtime = config["resources"]["metabat2"]["runtime"]
    shell:
        f"""
        mkdir -p {PROJECT_ROOT}/results/metabat2/{{wildcards.sample}}/depth
        jgi_summarize_bam_contig_depths --outputDepth {PROJECT_ROOT}/results/metabat2/{{wildcards.sample}}/depth/depth.txt {{input.bam}}
        metabat2 -i {{input.contigs}} -a {PROJECT_ROOT}/results/metabat2/{{wildcards.sample}}/depth/depth.txt \
                 -o {PROJECT_ROOT}/results/metabat2/{{wildcards.sample}}/bins/bin
        """

########################################
# 4a. CheckM: Quality control of bins + filtering
########################################
rule checkm:
    input:
        bins = lambda wc: os.path.join(PROJECT_ROOT, f"results/metabat2/{wc.sample}/bins")
    output:
        summary = os.path.join(PROJECT_ROOT, "results/checkm/{sample}/checkm_summary.tsv"),
        filtered_table = os.path.join(PROJECT_ROOT, "results/checkm/{sample}/checkm_filtered_bins.tsv"),
        filtered_bins = directory(os.path.join(PROJECT_ROOT, "results/checkm/{sample}/filtered_bins"))
    conda:
        lambda wc: os.path.join(PROJECT_ROOT, config["conda_envs"]["checkm"])
    threads: config["resources"]["checkm"]["threads"]
    resources:
        mem_mb = config["resources"]["checkm"]["mem_mb"],
        runtime = config["resources"]["checkm"]["runtime"]
    params:
        completion = config["resources"]["checkm"]["completion_threshold"],
        contamination = config["resources"]["checkm"]["contamination_threshold"]
    shell:
        f"""
        mkdir -p {PROJECT_ROOT}/results/checkm/{{wildcards.sample}}

        # Run CheckM lineage workflow
        checkm lineage_wf -x fa {{input.bins}} {PROJECT_ROOT}/results/checkm/{{wildcards.sample}} -t {{threads}}

        # Generate QA summary table
        checkm qa {PROJECT_ROOT}/results/checkm/{{wildcards.sample}}/lineage.ms \
                 {PROJECT_ROOT}/results/checkm/{{wildcards.sample}} \
                 -o 2 > {{output.summary}}

        # Filter bins by thresholds (keep header + passing MAGs)
        awk -v comp={{params.completion}} -v cont={{params.contamination}} \
            'NR==1 || ($12+0 >= comp && $13+0 <= cont)' {{output.summary}} \
            > {{output.filtered_table}}

        # Create directory for filtered bins
        mkdir -p {{output.filtered_bins}}

        # Extract passing bin names and copy them
        tail -n +2 {{output.filtered_table}} | cut -f1 | while read bin; do
            if [ -f {{input.bins}}/$bin.fa ]; then
                cp {{input.bins}}/$bin.fa {{output.filtered_bins}}/$bin.fa
            fi
        done
        """

########################################
# 4b. Summarize MAGs passing CheckM thresholds
########################################
rule summarize_checkm_hq:
    input:
        expand(os.path.join(PROJECT_ROOT, "results/checkm/{sample}/checkm_filtered_bins.tsv"), sample=SAMPLES)
    output:
        os.path.join(PROJECT_ROOT, "results/checkm/QC_summary.tsv")
    run:
        records = []
        for sample in SAMPLES:
            summary_path = os.path.join(PROJECT_ROOT, f"results/checkm/{sample}/checkm_summary.tsv")
            filtered_path = os.path.join(PROJECT_ROOT, f"results/checkm/{sample}/checkm_filtered_bins.tsv")
            if not os.path.exists(summary_path) or not os.path.exists(filtered_path):
                continue
            try:
                total_bins = sum(1 for line in open(summary_path)) - 1
                passed_bins = sum(1 for line in open(filtered_path)) - 1
                percent_passed = (passed_bins / total_bins * 100) if total_bins > 0 else 0
            except Exception:
                total_bins, passed_bins, percent_passed = 0, 0, 0
            records.append({
                "sample": sample,
                "total_bins": total_bins,
                "passed_bins": passed_bins,
                "percent_passed": round(percent_passed, 1)
            })
        df = pd.DataFrame(records)
        df.to_csv(output[0], sep="\t", index=False)

########################################
# 5. GTDB-Tk: Taxonomic classification (QC passing MAGs only)
########################################
rule gtdbtk:
    input:
        bins = lambda wc: os.path.join(PROJECT_ROOT, f"results/checkm/{wc.sample}/filtered_bins")
    output:
        classification = os.path.join(PROJECT_ROOT, "results/gtdbtk/{sample}/gtdbtk.bac120.summary.tsv")
    conda:
        os.path.join(PROJECT_ROOT, "envs/GTDBTK.yaml")
    threads: config["resources"]["gtdbtk"]["threads"]
    resources:
        mem_mb = config["resources"]["gtdbtk"]["mem_mb"],
        runtime = config["resources"]["gtdbtk"]["runtime"]
    params:
        gtdbtk_data = config["gtdbtk_data"]
    shell:
        f"""
        mkdir -p {PROJECT_ROOT}/results/gtdbtk/{{wildcards.sample}}
        gtdbtk classify_wf \
            --genome_dir {{input.bins}} \
            --out_dir {PROJECT_ROOT}/results/gtdbtk/{{wildcards.sample}} \
            --cpus {{threads}} \
            --data_dir {{params.gtdbtk_data}}
        """

########################################
# 5. CoverM: Abundance estimation (QC passing MAGs only)
########################################
rule coverm:
    input:
        bam = lambda wc: os.path.join(PROJECT_ROOT, f"results/samtools/{wc.sample}/{wc.sample}.sorted.bam"),
        bins = lambda wc: os.path.join(PROJECT_ROOT, f"results/checkm/{wc.sample}/filtered_bins")
    output:
        os.path.join(PROJECT_ROOT, "results/coverm/{sample}/abundance.tsv")
    conda:
        lambda wc: os.path.join(PROJECT_ROOT, config["conda_envs"]["coverm"])
    threads: config["resources"]["coverm"]["threads"]
    resources:
        mem_mb = config["resources"]["coverm"]["mem_mb"],
        runtime = config["resources"]["coverm"]["runtime"]
    shell:
        f"""
        mkdir -p {PROJECT_ROOT}/results/coverm/{{wildcards.sample}}
        coverm genome -b {{input.bam}} -m relative_abundance \
                      -s {{input.bins}} -t {{threads}} \
                      > {{output}}
        """
