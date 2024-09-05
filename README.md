# Bioinformatics Pipeline: Individual Scripts Version

## Overview

This pipeline consists of a series of individual scripts that perform various bioinformatics analyses, including read filtering, taxonomic classification, and BLAST searches. The pipeline is designed to process multiple samples in parallel, using a manifest file to specify input data.

## Components

1. `prepInput.sh`: Prepares input data and directory structure
2. `filterMerge.sh`: Filters and merges paired-end reads
3. `FqToFa.sh`: Converts FASTQ to FASTA format
4. `krakenUniq.sh`: Runs KrakenUniq for taxonomic classification
5. `sortKraken.R`: Processes and sorts Kraken results
6. `getSampKrakReads.sh`: Extracts classified reads from Kraken results
7. `splitReads.sh`: Splits reads into smaller files
8. `runBlastNt.sh`: Runs BLAST on the split reads
9. `blastFiltTopSamp.R`: Processes BLAST results
10. `tax.py`: Performs taxonomic analysis
11. `taxanalysis.py`: Calculates taxonomic metrics

## Prerequisites

- Bash
- R
- Python 3
- fastp
- seqtk
- KrakenUniq
- BLAST
- seqkit

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/your_username/bioinfo_pipeline.git
   cd bioinfo_pipeline
   ```

2. Install required software (see Prerequisites).

3. Set up necessary databases (KrakenUniq, BLAST nt).

## Usage

1. Prepare your manifest file (`manifest.tsv`) with three columns:
   - Path to read 1
   - Path to read 2
   - Sample prefix

2. Run each script in order, using the manifest file as input. For example:

   ```bash
   ./prepInput.sh manifest.tsv
   ./filterMerge.sh manifest.tsv /path/to/adapter.fasta
   ./FqToFa.sh manifest.tsv
   ./krakenUniq.sh manifest.tsv /path/to/krakenuniq_db
   # ... and so on for each script
   ```

3. The final results will be in the `output` directory.

## Output

- Filtered and merged reads
- Taxonomic classification results
- BLAST results
- Taxonomic analysis and metrics

## Notes

- Adjust the number of threads or parallel jobs in scripts based on your system's capabilities.
- Ensure you have sufficient disk space, especially for the BLAST database and intermediate files.

# Bioinformatics Pipeline: Combined Bash Script Version

## Overview

This version of the pipeline combines all individual scripts into a single Bash script, streamlining the execution process. The pipeline performs various bioinformatics analyses, including read filtering, taxonomic classification, and BLAST searches, processing multiple samples in parallel.

## Prerequisites

- Bash
- R
- Python 3
- fastp
- seqtk
- KrakenUniq
- BLAST
- seqkit

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/your_username/bioinfo_pipeline.git
   cd bioinfo_pipeline
   ```

2. Install required software (see Prerequisites).

3. Set up necessary databases (KrakenUniq, BLAST nt).

4. Make the script executable:
   ```
   chmod +x pipeline.sh
   ```

## Usage

Run the pipeline with the following command:

```bash
./pipeline.sh <manifest_file> <adapter_fasta> <krakenuniq_db> <blast_db> <output_dir>
```

Arguments:
- `<manifest_file>`: Path to the manifest file with columns: read1_path, read2_path, prefix
- `<adapter_fasta>`: Path to the adapter FASTA file for fastp
- `<krakenuniq_db>`: Path to the KrakenUniq database
- `<blast_db>`: Path to the BLAST database
- `<output_dir>`: Path to the output directory

Example:
```bash
./pipeline.sh manifest.tsv adapters.fa /path/to/kraken_db /path/to/blast_db results
```

## Pipeline Steps

1. Prepare input
2. Filter and merge reads
3. Convert FASTQ to FASTA
4. Run KrakenUniq
5. Sort Kraken results
6. Extract Kraken classified reads
7. Split reads
8. Run BLAST
9. Perform taxonomic analysis
10. Calculate taxonomic metrics

## Output

The script will create the specified output directory containing:
- Processed reads
- Taxonomic classification results
- BLAST results
- Taxonomic analysis files
- Pipeline performance metrics

## Notes

- Adjust the number of threads or parallel jobs in the script based on your system's capabilities.
- Ensure you have sufficient disk space, especially for the BLAST database and intermediate files.
- The script uses relative paths, so run it from the directory containing the script and input files.

  # Bioinformatics Pipeline: Nextflow Version

## Overview

This version of the pipeline is implemented using Nextflow, a powerful tool for building scalable and reproducible scientific workflows. The pipeline performs various bioinformatics analyses, including read filtering, taxonomic classification, and BLAST searches, with improved parallelization and portability.

## Prerequisites

- Nextflow
- Java 8 or later
- Docker or Singularity (optional, for containerization)

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/your_username/bioinfo_pipeline.git
   cd bioinfo_pipeline
   ```

2. Install Nextflow:
   ```
   curl -s https://get.nextflow.io | bash
   ```

3. (Optional) Pull the Docker image or build the Singularity image if using containers.

## Usage

Run the pipeline with the following command:

```bash
nextflow run main.nf --manifest path/to/manifest.tsv --adapter_fasta path/to/adapter.fasta --krakenuniq_db path/to/krakenuniq_db --blast_db path/to/blast_db --output_dir results
```

Parameters:
- `--manifest`: Path to the manifest file with columns: read1_path, read2_path, prefix
- `--adapter_fasta`: Path to the adapter FASTA file for fastp
- `--krakenuniq_db`: Path to the KrakenUniq database
- `--blast_db`: Path to the BLAST database
- `--output_dir`: Path to the output directory (default: results)

## Pipeline Processes

1. filterAndMerge: Filters and merges paired-end reads
2. fastqToFasta: Converts FASTQ to FASTA format
3. runKrakenUniq: Runs KrakenUniq for taxonomic classification
4. sortKraken: Processes and sorts Kraken results
5. extractKrakenReads: Extracts classified reads from Kraken results
6. splitReads: Splits reads into smaller files
7. runBlast: Runs BLAST on the split reads
8. taxonomicAnalysis: Performs detailed taxonomic analysis
9. calculateMetrics: Calculates precision metrics for the pipeline

## Configuration

You can customize the pipeline execution by modifying the `nextflow.config` file. This file allows you to specify:

- Executor (e.g., local, slurm, sge)
- Resource allocations (CPUs, memory, time)
- Container usage (Docker/Singularity)
- Pipeline parameters

Example `nextflow.config`:

```groovy
process {
    executor = 'slurm'
    queue = 'normal'
    cpus = 8
    memory = '16 GB'
    time = '2h'
}

params {
    manifest = "/path/to/your/manifest.tsv"
    adapter_fasta = "/path/to/your/adapter.fasta"
    krakenuniq_db = "/path/to/your/krakenuniq_db"
    blast_db = "/path/to/your/blast_db"
    output_dir = "results"
}
```

## Output

The pipeline will create the specified output directory containing:
- Processed reads
- Taxonomic classification results
- BLAST results
- Taxonomic analysis files
- Pipeline performance metrics

## Notes

- The Nextflow pipeline automatically handles parallelization and can be easily scaled to different computational environments.
- Use the `-resume` flag to restart the pipeline from the last successful step if it's interrupted.
- For large-scale analyses, consider using a container technology (Docker or Singularity) to ensure consistency across different environments.
