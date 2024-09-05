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
    echo "Extracting Kraken classified reads for $prefix"
    cd "$prefix"
    seqtk subseq "${prefix}_krakUniq_classified.reads" krakenReads.id > "${prefix}_selectKraken.reads"
    cd ..
done < "../$manifest_file"

cd ..
echo "Extraction of Kraken classified reads complete."
