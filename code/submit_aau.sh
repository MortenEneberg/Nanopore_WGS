#!/bin/sh
# Load all required modules for the job
module load snakemake/7.18.2-foss-2020b
module load Miniconda3/4.9.2-foss-2019a

WD="/user_data/men/sepseq/SepSeeQ/2023-10_WGS_Snakemake_pipeline/Nanopore_WGS/"                     #insert AAU WD.
cd $WD

#  Put your job commands after this line. Load all required modules before submitting this script.
snakemake --latency-wait 90 -s Snakefile --configfile config/config.yaml --cores 40 --use-conda --conda-frontend conda 
