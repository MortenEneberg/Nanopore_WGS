# Nanopore WGS

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.18.2-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow to assemble genomes from Nanopore WGS data

# Getting started
```
git clone https://github.com/MortenEneberg/Nanopore_WGS/
cd Nanopore_WGS
```

# Running the Snakemake pipeline

The pipeline is executed via the submission scripts (e.g. `bash code/submit_aau.sh`). Here the modules for loading conda and snakemake are specified. Alternatively, execute the pipeline in a conda env with these dependencies.

Specify paths to:
- gtdb-tk database in the gtdbtk.yaml environment file
- In the config file:
  - reads: This folder should be ordered such that it has library folders with subfolders of barcodes with .gz reads.
  - metadata: metadata connecting samples to seq library and barcode (a csv file with 3 columns: library, barcode, and sampleID)
  - checkm2db: checkm2 diamond database (https://github.com/chklovski/CheckM2#database)
  - reference_genome: optional
