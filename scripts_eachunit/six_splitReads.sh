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
    echo "Splitting reads for $prefix"
    cd "$prefix"
    rm -rf splitSeq10K
    mkdir splitSeq10K
    seqkit split "${prefix}_selectKraken.reads" -s 900 -j 8 -O ./splitSeq10K/
    cd splitSeq10K
    for j in *.reads; do mkdir "${j}_dir"; mv "$j" "./${j}_dir/reads.fas"; done
    cd ../..
done < "../$manifest_file"

cd ..
echo "Read splitting complete."
