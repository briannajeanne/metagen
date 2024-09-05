#!/bin/bash

# Check if manifest file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <path_to_manifest_file>"
    exit 1
fi

manifest_file=$1

# Create process directory if it doesn't exist
mkdir -p process

# Read the manifest file and process each line
while IFS=$'\t' read -r read1 read2 prefix
do
    # Create directory for the sample
    mkdir -p "./process/$prefix"
    
    # Copy input files to the sample directory
    cp "$read1" "./process/$prefix/${prefix}_R1.fastq.gz"
    cp "$read2" "./process/$prefix/${prefix}_R2.fastq.gz"
    
    # Create sample.name file
    echo "$prefix" > "./process/$prefix/sample.name"
    
    # Add to newdir.list
    echo "$prefix" >> newdir.list
done < "$manifest_file"

echo "Input preparation complete."
