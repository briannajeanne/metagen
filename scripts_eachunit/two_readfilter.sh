#!/bin/bash

# Check if manifest file is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <path_to_manifest_file> <adapter_fasta_path>"
    exit 1
fi

manifest_file=$1
adapter_fasta=$2

cd process

while IFS=$'\t' read -r read1 read2 prefix
do
    echo "Processing $prefix"
    cd "$prefix"
    fastp -i "${prefix}_R1.fastq.gz" -I "${prefix}_R2.fastq.gz" \
    --detect_adapter_for_pe \
    --adapter_fasta "$adapter_fasta" \
    --merge --merged_out "${prefix}_merged_reads.fastq" \
    --include_unmerged \
    -l 100 \
    --thread 8
    cd ..
done < "../$manifest_file"

cd ..
echo "Filtering and merging complete."
