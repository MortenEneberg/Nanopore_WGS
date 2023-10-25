import string
import sys
import pathlib
import os
import subprocess
import gzip
import pwd
import datetime
import pandas as pd
import shutil
import csv
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

def process_sample_fastq_file(metadata_path_set, root_folder_path_set, output_dir_set, target_sample_id_set):
	"""
    Process .fastq.gz files based on the provided sample ID, extract information from metadata,
    and handle the files accordingly.

    :param metadata_path: Path to the metadata CSV file.
    :param read_folder_path: Base path of the read folders.
    :param output_folder: Destination to save the processed files.
    :param sample_id: The unique sample ID to process.
	"""
	metadata_path = next(iter(metadata_path_set))
	root_folder_path = next(iter(root_folder_path_set))
	output_dir = next(iter(output_dir_set))
	sample_id = next(iter(target_sample_id_set))

	sample_id_found = False  # Flag to check if the sample was found in metadata.
	processed_files = []  # Keeping track of which files have been processed.

	try:
		with open(metadata_path, 'r', encoding='utf-8-sig') as infile:
			reader = csv.DictReader(infile)

            # Searching for the sample ID in the metadata.	
			for row in reader:
				current_sample_id = row.get('sampleID')
				if current_sample_id == sample_id:
					sample_id_found = True
					lib = row['lib']
					bar = row['bar']

					# Construct the file path for .fastq.gz files.
					fastq_file_path = os.path.join(root_folder_path, lib, bar)

					# Check if the directory exists before proceeding.
					if not os.path.exists(fastq_file_path):
						print(f"Directory does not exist: {fastq_file_path}")
						continue

					 # Prepare the output directory.
					sample_output_folder = os.path.join(output_dir, sample_id)
					os.makedirs(sample_output_folder, exist_ok=True)

					# For storing the content of .fastq files.
					concatenated_content = []

					# Process the .fastq.gz files.
					for file_name in os.listdir(fastq_file_path):
						if file_name.endswith('.fastq.gz'):
							# Full path to the .fastq.gz file.
							full_file_path = os.path.join(fastq_file_path, file_name)
							
							# Read the contents of the .fastq.gz file.
							with gzip.open(full_file_path, 'rt') as f:  # 'rt' mode for reading text.
								content = f.read()
								concatenated_content.append(content)

							processed_files.append(full_file_path)

					# Concatenate contents, write to a new .fastq file, and compress it.
					if concatenated_content:
						output_fastq_file = os.path.join(sample_output_folder, f"{sample_id}_reads.fastq")
						with open(output_fastq_file, 'w') as f_out:
							f_out.write(''.join(concatenated_content))

					break

			if not sample_id_found:
				print(f"Sample ID {sample_id} not found in the metadata.")
			elif processed_files:
				print(f"Processed files for sample ID {sample_id}: {', '.join(processed_files)}")
			else:
				print(f"No .fastq.gz files found for sample ID {sample_id} in the expected directory.")

	except Exception as e:
		print(f"An error occurred: {e}")

    # You can return something from your function if needed, for example, the path to the new file.
	return os.path.join(output_dir, sample_id, f"{sample_id}_reads.fastq.gz") if sample_id_found else None


def extract_sampleID_column(file_path):
    # List to store the extracted data
    sample_ids = []

    # Open the CSV file
    with open(file_path, mode='r', newline='', encoding='utf-8') as csvfile:
        # Use the CSV DictReader to easily access columns by name
        reader = csv.DictReader(csvfile)

        # Check if 'sampleID' exists in the file
        if 'sampleID' in reader.fieldnames:
            # Iterate over rows and add the 'sampleID' values to our list
            for row in reader:
                sample_ids.append(row['sampleID'])
        else:
            raise ValueError("No 'sampleID' column in the CSV")

    return sample_ids

READ_FOLDER = config["reads"]
METADATA = config["metadata"]
sample_ids = extract_sampleID_column(METADATA)
minQual = config["minReadsQual"]
minLength = config["minReadsLength"]
 
##### Define Rules !

rule all:
	input:
		expand("data/20_flye_assembly/{sampleid}/assembly.fasta",sampleid=sample_ids),
		expand("data/qc/01_NanoPlot_raw/{sampleid}/NanoPlot-report.html",sampleid=sample_ids),
		expand("data/qc/10_quast/{sampleid}/report.html",sampleid=sample_ids),
		expand("data/qc/20_checkm/{sampleid}/{sampleid}_CheckM.txt",sampleid=sample_ids),
		expand("data/qc/02_NanoPlot_filtered/{sampleid}/NanoPlot-report.html",sampleid=sample_ids),
		expand("data/AMR/01_abricate/{sampleid}/abricate.txt",sampleid=sample_ids)
		#expand("/{smp}/ResistanceAbricate/{smp}_abricate.txt",smp=SAMPLES),
		#expand(outputdir+"/{smp}/BuscoOutput/{smp}_busco.txt",smp=SAMPLES),
		#expand(outputdir+"/{smp}/BuscoOutput/{smp}_busco_parsed.txt",smp=SAMPLES),
		#expand(outputdir+"/{smp}/Mob_recon/contig_report.txt",smp=SAMPLES),
		#expand(outputdir+"/{smp}/Prokka/{smp}.gff",smp=SAMPLES),
		#outputdir+"/roary/core_gene_alignment.aln",
		#outputdir+"/FastTree/tree.newick"

rule cat_raw_reads:
	output:
		fastq = "data/01_cat_reads/{sampleid}/{sampleid}_reads.fastq"
	params:
		outdir = "data/01_cat_reads/",
		reads = config["reads"],
		metadata = config["metadata"],
	run:
		process_sample_fastq_file({params.metadata}, {params.reads}, {params.outdir}, {wildcards.sampleid})

rule Nanoplot_raw:
	conda:
		"_envs/nanoplot.yaml"
	input: 
		"data/01_cat_reads/{sampleid}/{sampleid}_reads.fastq"
	output:
		"data/qc/01_NanoPlot_raw/{sampleid}/NanoPlot-report.html"
	params:
		outdir_nanoplot = "data/qc/01_NanoPlot_raw"
	shell:
		"""
		NanoPlot --plots dot -f svg --N50 --fastq {input} -o {params.outdir_nanoplot}/{wildcards.sampleid} 
		"""

rule filter_reads:
	conda:
		"_envs/filtlong.yaml"
	input:
		"data/01_cat_reads/{sampleid}/{sampleid}_reads.fastq"
	output:
		"data/10_filtered_reads/{sampleid}/{sampleid}_filtered.fastq"
	params:
		minl = minLength,
		minqual = minQual
	shell:
		"""
		filtlong --min_length {params.minl} --min_mean_q {params.minqual} {input} > {output}
		"""

rule Nanoplot_filtered:
	conda:
		"_envs/nanoplot.yaml"
	input: 
		"data/10_filtered_reads/{sampleid}/{sampleid}_filtered.fastq"
	output:
		"data/qc/02_NanoPlot_filtered/{sampleid}/NanoPlot-report.html"
	params:
		outdir_nanoplot = "data/qc/02_NanoPlot_filtered"
	shell:
		"""
		NanoPlot --plots dot -f svg --N50 --fastq {input} -o {params.outdir_nanoplot}/{wildcards.sampleid} 
		"""

rule flye:
	conda:
		"_envs/flye.yaml"
	input:
		"data/10_filtered_reads/{sampleid}/{sampleid}_filtered.fastq"
	output:
		"data/20_flye_assembly/{sampleid}/assembly.fasta"
	params:
		outdir_flye = "data/20_flye_assembly",
		threads = config["threads_Flye"]
	shell:
		"""
		flye --nano-hq {input} --threads {params.threads} --out-dir {params.outdir_flye}/{wildcards.sampleid}
		"""

rule quast:
	conda:
		"_envs/quast.yaml"
	input:
		assembly = "data/20_flye_assembly/{sampleid}/assembly.fasta",
		reads = "data/10_filtered_reads/{sampleid}/{sampleid}_filtered.fastq"
	output:
		"data/qc/10_quast/{sampleid}/report.html"
	params:
		outdir = "data/qc/10_quast/",
		ref = config["reference_genome"]
	threads: 
		5
	shell:
		"""
		quast --threads {threads} --conserved-genes-finding --output-dir {params.outdir}/{wildcards.sampleid} -r {params.ref} -L {input.assembly} --nanopore {input.reads}
        """

rule checkm:
	conda:
		"_envs/checkm.yaml"
	input:
		"data/20_flye_assembly/{sampleid}/assembly.fasta"
	output:
		"data/qc/20_checkm/{sampleid}/{sampleid}_CheckM.txt"
	params:
		assembly_folder = "data/20_flye_assembly/",
		outdir = "data/qc/20_checkm"
	threads:
		10
	shell:
		"""
		checkm lineage_wf -f {output} -t {threads} -x fasta {params.assembly_folder}/{wildcards.sampleid} {params.outdir}/{wildcards.sampleid}/
		"""

rule abricate:
	conda:
		"_envs/abricate.yaml"
	input:
		"data/20_flye_assembly/{sampleid}/assembly.fasta"
	output:
		"data/AMR/01_abricate/{sampleid}/abricate.txt"
	shell:
		"""
		abricate {input} > {output}
		"""