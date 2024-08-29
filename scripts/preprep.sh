#!/bin/bash

# Check if required arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <path_to_R1> <path_to_R2>"
    exit 1
fi

# Assign input arguments to variables
R1_PATH="$1"
R2_PATH="$2"

# Create input directory if it doesn't exist
mkdir -p input

# Loop through each fasta file in the current directory and generate reads using iss
for fasta in *.fasta; do
    # Check if there are any .fasta files
    if [ ! -f "$fasta" ]; then
        echo "No .fasta files found in the current directory."
        exit 1
    fi

    # Get the base name of the fasta file (without extension)
    base_name=$(basename "$fasta" .fasta)
    
    # Generate reads with iss
    iss generate -g "$fasta" -o "input/${base_name}_reads" -n 3000 --cpus 16 --model miseq
done

# Check if the R script exists before running it
if [ -f "./programs/bin/create_library.R" ]; then
    # Create the library CSV using the R script
    Rscript --vanilla ./programs/bin/create_library.R
else
    echo "Warning: create_library.R script not found. Skipping library creation."
fi

# Concatenate all _R1 reads into one combined file
cat input/*_R1.fastq > "$R1_PATH"

# Concatenate all _R2 reads into one combined file
cat input/*_R2.fastq > "$R2_PATH"

# Gzip both files
gzip "$R1_PATH"
gzip "$R2_PATH"

echo "Processing complete. Output files: ${R1_PATH}.gz and ${R2_PATH}.gz"
