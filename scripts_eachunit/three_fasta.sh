#!/bin/bash

# Check if manifest file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <path_to_manifest_file>"
    exit 1
fi

manifest_file=$1

cd process

while IFS=$'\t' read -r read1 read2 prefix
do
    echo "Converting $prefix to FASTA"
    cd "$prefix"
    seqtk seq -a "${prefix}_merged_reads.fastq" > "${prefix}_merged_reads.fa"
    cd ..
done < "../$manifest_file"

cd ..
echo "FASTQ to FASTA conversion complete."
