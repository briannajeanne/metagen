#!/bin/bash

# Check if manifest file and database path are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <path_to_manifest_file> <path_to_krakenuniq_db>"
    exit 1
fi

manifest_file=$1
krakenuniq_db=$2

cd process

while IFS=$'\t' read -r read1 read2 prefix
do
    echo "Running KrakenUniq on $prefix"
    cd "$prefix"
    krakenuniq --threads 8 --db "$krakenuniq_db" \
    --report-file "${prefix}_krakUniq_sample.report" --output "${prefix}_krakUniq_sample.kraken" \
    --classified-out "${prefix}_krakUniq_classified.reads" "${prefix}_merged_reads.fa"
    cd ..
done < "../$manifest_file"

cd ..
echo "KrakenUniq classification complete."
