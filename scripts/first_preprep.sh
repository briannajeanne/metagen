#!/bin/bash

# Check if the manifest file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <path_to_manifest_file>"
    exit 1
fi

MANIFEST_FILE="$1"

# Check if the manifest file exists
if [ ! -f "$MANIFEST_FILE" ]; then
    echo "Error: Manifest file not found: $MANIFEST_FILE"
    exit 1
fi

# Create input and output directories if they don't exist
mkdir -p input output

# Process each line in the manifest file
while IFS=$'\t' read -r R1_PATH R2_PATH || [ -n "$R1_PATH" ]; do
    # Skip empty lines
    [ -z "$R1_PATH" ] && continue

    # Check if both R1 and R2 paths are provided
    if [ -z "$R2_PATH" ]; then
        echo "Error: R2 path is missing for R1: $R1_PATH"
        continue
    fi

    # Check if the input files exist
    if [ ! -f "$R1_PATH" ] || [ ! -f "$R2_PATH" ]; then
        echo "Error: Input file(s) not found for R1: $R1_PATH or R2: $R2_PATH"
        continue
    fi

    # Get the base name for the output
    BASE_NAME=$(basename "$R1_PATH" _R1.fastq)

    # Copy input files to the input directory
    cp "$R1_PATH" "input/${BASE_NAME}_R1.fastq"
    cp "$R2_PATH" "input/${BASE_NAME}_R2.fastq"

    # Generate reads with iss (assuming we're using the input FASTQ files as templates)
    iss generate --genomes "input/${BASE_NAME}_R1.fastq" "input/${BASE_NAME}_R2.fastq" \
                 -o "output/${BASE_NAME}_reads" -n 3000 --cpus 16 --model miseq

    echo "Processed: $BASE_NAME"

done < "$MANIFEST_FILE"

# Check if the R script exists before running it
if [ -f "./programs/bin/create_library.R" ]; then
    # Create the library CSV using the R script
    Rscript --vanilla ./programs/bin/create_library.R
else
    echo "Warning: create_library.R script not found. Skipping library creation."
fi

# Concatenate all generated R1 reads into one combined file
cat output/*_R1.fastq > output/R1_001.fastq

# Concatenate all generated R2 reads into one combined file
cat output/*_R2.fastq > output/R2_001.fastq

# Gzip both files
gzip output/R1_001.fastq
gzip output/R2_001.fastq

echo "Processing complete. Output files: output/R1_001.fastq.gz and output/R2_001.fastq.gz"
