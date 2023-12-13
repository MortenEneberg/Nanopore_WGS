#!/usr/bin/env python

import os
import argparse
from random import sample
import subprocess
import string
import sys
import pathlib
import gzip
import pwd
import datetime
import shutil
import csv
import re  # You were missing the import for the 're' module

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--metadata", type=str)
parser.add_argument("-o", "--outdir", type=str)
parser.add_argument("-r", "--reads", type=str)
parser.add_argument("-s", "--sample", type=str)

args = parser.parse_args()
metadata = args.metadata
outdir = args.outdir
reads = args.reads
sample = args.sample



def process_sample_fastq_file(metadata_path_set, read_folder_path_set, output_dir_set, target_sample_id_set):
	"""
    Process .fastq.gz files based on the provided sample ID, extract information from metadata,
    and handle the files accordingly.

    :param metadata_path: Path to the metadata CSV file.
    :param read_folder_path: Base path of the read folders.
    :param output_folder: Destination to save the processed files.
    :param sample_id: The unique sample ID to process.
	"""
	metadata_path = next(iter(metadata_path_set))
	read_folder_path = next(iter(read_folder_path_set))
	output_folder = next(iter(output_dir_set))
	sample_id = next(iter(target_sample_id_set))

	sample_id_found = False  # Flag to check if the sample was found in metadata.
	processed_files = []  # Keeping track of which files have been processed.

	try:
		with open(metadata_path, 'r', encoding='utf-8-sig') as infile:
			reader = csv.DictReader(infile)

			for row in reader:
				# Adjusted column names based on your updates
				if row.get('sample_id') == sample_id:
					sample_id_found = True
					lib_flowcell_id = row['lib_flowcell_id']
					lib_barcode = row['lib_barcode']
					print(lib_barcode)
					print(lib_flowcell_id)
					# Construct the file path for .fastq.gz files.
					fastq_file_path = os.path.join(read_folder_path, lib_flowcell_id, lib_barcode)

					# Check if the directory exists before proceeding.
					if not os.path.exists(fastq_file_path):
						print(f"Directory does not exist: {fastq_file_path}")
						continue

					# Prepare the output directory.
					sample_output_folder = os.path.join(output_folder, sample_id)
					os.makedirs(sample_output_folder, exist_ok=True)

					# For storing the content of .fastq files.
					concatenated_content = []

					# Process the .fastq.gz files.
					for file_name in os.listdir(fastq_file_path):
						if file_name.endswith('.fastq.gz'):
							full_file_path = os.path.join(fastq_file_path, file_name)
							
							with gzip.open(full_file_path, 'rt') as f:  # 'rt' mode for reading text.
								content = f.read()
								concatenated_content.append(content)

							processed_files.append(full_file_path)

					# Concatenate contents, write to a new .fastq file, and compress it.
					if concatenated_content:
						output_fastq_file = os.path.join(sample_output_folder, f"{sample_id}.fastq")
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

	return os.path.join(output_folder, sample_id, f"{sample_id}.fastq.gz") if sample_id_found else None

process_sample_fastq_file(metadata, reads, outdir, sample)