#!/bin/bash

# Check if required arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <manifest_file> <adapter_fasta_path>"
    exit 1
fi

MANIFEST_FILE="$1"
ADAPTER_FASTA="$2"

# Check if manifest file exists
if [ ! -f "$MANIFEST_FILE" ]; then
    echo "Error: Manifest file not found: $MANIFEST_FILE"
    exit 1
fi

# Check if adapter FASTA file exists
if [ ! -f "$ADAPTER_FASTA" ]; then
    echo "Error: Adapter FASTA file not found: $ADAPTER_FASTA"
    exit 1
fi

# Create a process directory if it doesn't exist
PROCESS_DIR="process"
mkdir -p "$PROCESS_DIR"

# Move to the process directory
cd "$PROCESS_DIR" || exit 1

# Process each line in the manifest file
while IFS=$'\t' read -r R1_PATH R2_PATH OUTPUT_PREFIX || [ -n "$R1_PATH" ]; do
    # Skip empty lines
    [ -z "$R1_PATH" ] && continue

    # Check if all three columns are provided
    if [ -z "$R2_PATH" ] || [ -z "$OUTPUT_PREFIX" ]; then
        echo "Error: Missing R2 path or output prefix for R1: $R1_PATH"
        continue
    fi

    # Check if the input files exist
    if [ ! -f "$R1_PATH" ] || [ ! -f "$R2_PATH" ]; then
        echo "Error: Input file(s) not found for R1: $R1_PATH or R2: $R2_PATH"
        continue
    fi

    # Run fastp on the paired-end files
    fastp \
        -i "$R1_PATH" \
        -I "$R2_PATH" \
        --detect_adapter_for_pe \
        --adapter_fasta "$ADAPTER_FASTA" \
        --merge \
        --merged_out "${OUTPUT_PREFIX}_merged_reads.fastq" \
        --include_unmerged \
        -l 100 \
        --thread 8 \
        -o "${OUTPUT_PREFIX}_R1_trimmed.fastq.gz" \
        -O "${OUTPUT_PREFIX}_R2_trimmed.fastq.gz" \
        -j "${OUTPUT_PREFIX}_fastp.json" \
        -h "${OUTPUT_PREFIX}_fastp.html"

    echo "Processed: $OUTPUT_PREFIX"

done < "../$MANIFEST_FILE"

echo "FastP processing complete. Output files are in the '$PROCESS_DIR' directory."
