#!/bin/bash

# Check if manifest file and BLAST database path are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <path_to_manifest_file> <path_to_blast_db>"
    exit 1
fi

manifest_file=$1
blast_db=$2

while IFS=$'\t' read -r read1 read2 prefix
do
    echo "Running BLAST for $prefix"
    cd "process/$prefix/splitSeq10K"
    ls -d */ | parallel -j 28 "cd {} && blastn -db $blast_db -query reads.fas \
    -num_threads 3 -out ${prefix}_NtV4_blast.results \
    -outfmt '7 staxids qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 10"
    cd ../../../
done < "$manifest_file"

echo "BLAST searches complete."
