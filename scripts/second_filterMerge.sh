#!/bin/bash

# Check if required arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <adapter_fasta_path>"
    exit 1
fi

INPUT_DIR="$1"
ADAPTER_FASTA="$2"

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory not found: $INPUT_DIR"
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

# Create symbolic links to input files
ln -s "../$INPUT_DIR"/*.fastq.gz .

# Run fastp on all paired-end files
ls *R1_001.fastq.gz | sed 's/_R1_001.fastq.gz$//' | parallel -j 4 'fastp \
    -i {}_R1_001.fastq.gz \
    -I {}_R2_001.fastq.gz \
    --detect_adapter_for_pe \
    --adapter_fasta '"$ADAPTER_FASTA"' \
    --merge \
    --merged_out {}_merged_reads.fastq \
    --include_unmerged \
    -l 100 \
    --thread 8 \
    -o {}_R1_trimmed.fastq.gz \
    -O {}_R2_trimmed.fastq.gz \
    -j {}_fastp.json \
    -h {}_fastp.html'

echo "FastP processing complete. Output files are in the '$PROCESS_DIR' directory."
