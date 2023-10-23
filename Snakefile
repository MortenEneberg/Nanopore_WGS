import string
import sys
import pathlib
import os
import subprocess
import gzip
import pwd
import datetime
import re  # You were missing the import for the 're' module

def getCurrentTime():
    return datetime.datetime.now().strftime("%Y_%b_%d_%H:%M")

def getCurrentDate():
    return datetime.datetime.now().strftime("%y%m%d")

def GetFilename(FullPath):
    return os.path.splitext(os.path.basename(FullPath))[0]

def FixDirectoryName(Dir):
    x = re.search('/$', Dir)
    if x is not None:
        return Dir
    else:
        return Dir + "/"

def unique(list1):
    unique_list = []
    for x in list1:
        if x not in unique_list:
            unique_list.append(x)
    return unique_list

def generate_barcodes():
    barcodes = ["barcode" + str(i) for i in range(1, 97)]
    return barcodes

def list_directories(folder_path):
    directories = [d for d in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, d))]
    return directories

READ_FOLDER = config["reads"]
LIBRARIES = list_directories(READ_FOLDER)
BARCODES = generate_barcodes()
minQual = config["minReadsQual"]
minLength = config["minReadsLength"]
 
##### Define Rules !

rule all:
	input:
		expand("data/30_polished_assembly/{seq_lib}/{barcode}/assembly_{barcode}.fasta",seq_lib=LIBRARIES, barcode=BARCODES),
		expand("data/qc/01_NanoPlot_raw/{seq_lib}/{barcode}/{seq_lib}_{barcode}.html",seq_lib=LIBRARIES, barcode=BARCODES),
		#expand("/{smp}/ResistanceAbricate/{smp}_abricate.txt",smp=SAMPLES),
		#expand(outputdir+"/{smp}/BuscoOutput/{smp}_busco.txt",smp=SAMPLES),
		#expand(outputdir+"/{smp}/BuscoOutput/{smp}_busco_parsed.txt",smp=SAMPLES),
		#expand(outputdir+"/{smp}/Mob_recon/contig_report.txt",smp=SAMPLES),
		#expand(outputdir+"/{smp}/Prokka/{smp}.gff",smp=SAMPLES),
		#outputdir+"/roary/core_gene_alignment.aln",
		#outputdir+"/FastTree/tree.newick"

rule cat_raw_reads:
    output:
        fastq = "data/01_cat_reads/{seq_lib}/{barcode}.fastq"
    params:
        reads = READ_FOLDER
    shell:
        """
        cat {params.reads}{wildcards.seq_lib}/{wildcards.barcode}*.fastq.gz | gunzip -c > {output.fastq}
        """

rule Nanoplot_raw:
	conda:
		"_envs/nanoplot.yaml"
	input: 
		"data/01_cat_reads/{seq_lib}/{barcode}.fastq"
	output:
		"data/qc/01_NanoPlot_raw/{seq_lib}/{barcode}/{seq_lib}_{barcode}.html"
	params:
		outdir_nanoplot = "data/qc/01_NanoPlot_raw/{seq_lib}/{barcode}"
	shell:
		"""
		NanoPlot --plots dot --title {wildcards.seq_lib}_{wildcards.barcode} -f svg --N50 --fastq {input} -o {params.outdir_nanoplot} 
		"""

rule filter_reads:
	conda:
		"_envs/filtlong.yaml"
	input:
		"data/01_cat_reads/{seq_lib}/{barcode}.fastq"
	output:
		"data/10_filtered_reads/{seq_lib}/{barcode}.fastq"
	params:
		minl = minLength,
		minqual = minQual
	shell:
		"""
		filtlong --min_length {params.minl} --min_mean_q {params.minqual} {input} > {output}
		"""

rule flye:
	conda:
		"_envs/flye.yaml"
	input:
		"data/10_filtered_reads/{seq_lib}/{barcode}.fastq" 
	output:
		"data/20_flye_assembly/{seq_lib}/{barcode}/assembly_{barcode}.fasta"
	params:
		outdir_flye = "data/20_flye_assembly/{seq_lib}/{barcode}/",
		threads = config["threads_Flye"]
	shell:
		"""
		flye --nano-hq {input} --threads {params.threads} --out-dir {params.outdir_flye}
		"""

rule polish:
	conda:
		"_envs/medaka.yaml"
	input:
		reads = "data/10_filtered_reads/{seq_lib}/{barcode}.fastq",
		assembly = "data/20_flye_assembly/{seq_lib}/{barcode}/assembly_{barcode}.fasta"
	output:
		"data/30_polished_assembly/{seq_lib}/{barcode}/assembly_{barcode}.fasta"
	params:
		medaCPU = config["threads_Medaka"],
		medaMod = config["medaka_model"],
		medaBatch = config["batch_Medaka"],
		outdir_medaka = "data/30_polished_assembly/{seq_lib}/{barcode}/",
	shell:
		"""
		medaka_consensus -i {input.reads}  -d {input.assembly} -t {params.medaCPU} -f -m {params.medaMod} -b {params.medaBatch}
		"""		


rule abricate:
	conda:
		"_envs/abricate.yaml"
	input:
		"data/30_polished_assembly/{seq_lib}/{barcode}/assembly_{barcode}.fasta"
	output:
		"data/AMR/01_abricate/{seq_lib}/{barcode}/{barcode}_abricate.txt"
	shell:
		"""
		abricate {input} > {output}
		"""
