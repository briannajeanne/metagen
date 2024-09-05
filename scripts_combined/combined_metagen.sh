#!/bin/bash

# Usage function
usage() {
    echo "Usage: $0 <manifest_file> <adapter_fasta> <krakenuniq_db> <blast_db> <output_dir>"
    echo "  <manifest_file>: Path to the manifest file with columns: read1_path, read2_path, prefix"
    echo "  <adapter_fasta>: Path to the adapter FASTA file for fastp"
    echo "  <krakenuniq_db>: Path to the KrakenUniq database"
    echo "  <blast_db>: Path to the BLAST database"
    echo "  <output_dir>: Path to the output directory"
    exit 1
}

# Check if correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    usage
fi

# Assign arguments to variables
manifest_file=$1
adapter_fasta=$2
krakenuniq_db=$3
blast_db=$4
output_dir=$5

# Create output directory
mkdir -p "$output_dir"

# Function to run R script
run_r_script() {
    Rscript - <<EOF
$1
EOF
}

# Function to run Python script
run_python_script() {
    python - <<EOF
$1
EOF
}

# 1. Prepare input
echo "Step 1: Preparing input"
mkdir -p process
while IFS=$'\t' read -r read1 read2 prefix
do
    mkdir -p "./process/$prefix"
    cp "$read1" "./process/$prefix/${prefix}_R1.fastq.gz"
    cp "$read2" "./process/$prefix/${prefix}_R2.fastq.gz"
    echo "$prefix" > "./process/$prefix/sample.name"
    echo "$prefix" >> newdir.list
done < "$manifest_file"

# 2. Filter and merge reads
echo "Step 2: Filtering and merging reads"
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

# 3. Convert FASTQ to FASTA
echo "Step 3: Converting FASTQ to FASTA"
cd process
while IFS=$'\t' read -r read1 read2 prefix
do
    echo "Converting $prefix to FASTA"
    cd "$prefix"
    seqtk seq -a "${prefix}_merged_reads.fastq" > "${prefix}_merged_reads.fa"
    cd ..
done < "../$manifest_file"
cd ..

# 4. Run KrakenUniq
echo "Step 4: Running KrakenUniq"
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

# 5. Sort Kraken results
echo "Step 5: Sorting Kraken results"
run_r_script "
library(tidyverse)

manifest <- read.delim('$manifest_file', header = FALSE, col.names = c('read1', 'read2', 'prefix'))

for (i in 1:nrow(manifest)) {
    prefix <- manifest\$prefix[i]
    setwd(file.path('process', prefix))
    
    sampleName <- prefix
    
    virusId <- function() {
        df <- read.delim(paste0(prefix, '_krakUniq_sample.report'), comment.char = '#')
        virus <- as.data.frame(df)
        colnames(virus)[colnames(virus) == 'taxName'] <- 'Virus'
        virusID <- data.frame(unique(virus[, c('taxID', 'Virus')]))
        colnames(virusID)[1] <- 'krakTax'
        return(virusID)
    }
    
    idList <- virusId()
    readNameTab <- read.delim(paste0(prefix, '_krakUniq_sample.kraken'), F, sep = '\t')
    colnames(readNameTab)[3] <- 'krakTax'
    colnames(readNameTab)[2] <- 'Read'
    selCol <- readNameTab[, c('krakTax', 'Read')]
    cReads <- plyr::join(idList, selCol, by = 'krakTax', type = 'left', match = 'all')
    allKrakenReads <- unique(cReads[, c('Virus', 'Read')])
    allKrakenReads\$Sample <- sampleName
    
    write.table(allKrakenReads, paste0(prefix, '_krakenSelVirReads.tsv'), row.names = F, sep = '\t', quote = F)
    
    reads <- unique(allKrakenReads[, 'Read'])
    write.table(reads, paste0(prefix, '_krakenReads.id'), row.names = F, col.names = F, quote = F)
    
    file_path <- paste0(prefix, '_krakUniq_sample.kraken')
    data <- read_lines(file_path)
    
    parse_line <- function(line) {
        elements <- strsplit(line, '\\s+')[[1]]
        Read <- elements[2]
        primary_tax_id <- as.numeric(elements[3])
        tax_ids_all <- str_extract_all(line, '\\s[0-9]+(?=:|$)')[[1]]
        tax_ids_all <- as.numeric(gsub('\\s+', '', tax_ids_all))
        if (primary_tax_id <= 1000000000 && primary_tax_id > 0 && !primary_tax_id %in% tax_ids_all) {
            tax_ids_all <- c(primary_tax_id, tax_ids_all)
        }
        valid_tax_ids <- unique(tax_ids_all[tax_ids_all <= 1000000000 & tax_ids_all > 0])
        top_tax_ids <- head(valid_tax_ids, 3)
        if (length(top_tax_ids) < 3) {
            top_tax_ids <- c(top_tax_ids, rep(NA, 3 - length(top_tax_ids)))
        }
        return(data.frame(Read = Read, KID1 = top_tax_ids[1], KID2 = top_tax_ids[2], KID3 = top_tax_ids[3]))
    }
    
    result <- bind_rows(lapply(data, parse_line))
    write_csv(result, paste0(prefix, '_krakenSelVirReads.tsv'))
    
    setwd('../..')
}
"

# 6. Get sample Kraken reads
echo "Step 6: Extracting Kraken classified reads"
cd process
while IFS=$'\t' read -r read1 read2 prefix
do
    echo "Extracting Kraken classified reads for $prefix"
    cd "$prefix"
    seqtk subseq "${prefix}_krakUniq_classified.reads" "${prefix}_krakenReads.id" > "${prefix}_selectKraken.reads"
    cd ..
done < "../$manifest_file"
cd ..

# 7. Split reads
echo "Step 7: Splitting reads"
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

# 8. Run BLAST
echo "Step 8: Running BLAST"
while IFS=$'\t' read -r read1 read2 prefix
do
    echo "Running BLAST for $prefix"
    cd "process/$prefix/splitSeq10K"
    ls -d */ | parallel -j 28 "cd {} && blastn -db $blast_db -query reads.fas \
    -num_threads 3 -out ${prefix}_NtV4_blast.results \
    -outfmt '7 staxids qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 10"
    cd ../../../
done < "$manifest_file"

# 9. Taxonomic analysis
echo "Step 9: Performing taxonomic analysis"
run_python_script "
import os
import sys
from ete3 import NCBITaxa
import pandas as pd
import itertools
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict

manifest_file = '$manifest_file'
output_dir = '$output_dir'

# Read manifest file
manifest = pd.read_csv(manifest_file, sep='\t', header=None, names=['read1', 'read2', 'prefix'])

# Initialize NCBITaxa
ncbi = NCBITaxa()

def get_rank(taxid):
    return ncbi.get_rank([taxid]).get(taxid, None)

def is_classified(taxid):
    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    return 'unclassified' not in ' '.join(names.values()).lower()

def get_lca(taxids):
    try:
        combinations = list(itertools.combinations(taxids, 2))
        lcas = [(pair, ncbi.get_topology(pair)) for pair in combinations if pair[0] and pair[1]]
        if not lcas:
            return None, None, None
        lcas.sort(key=lambda x: ncbi.get_rank([x[1].taxid])[x[1].taxid])
        return lcas[0][0][0], lcas[0][0][1], lcas[0][1].taxid
    except KeyError as e:
        print(f'KeyError encountered: {e}')
        return None, None, None

def process_row(row):
    bid_used, kid_used, lca_taxid = determine_final_lca(row)
    return row['Main_Read'], bid_used, kid_used, lca_taxid

def preprocess_dataframe(df):
    df['Main_Read'] = df['Read'].str.extract(r'(.*?)(?:/\d+)?$')
    df_1 = df[df['Read'].str.endswith('/1')].add_suffix('_1')
    df_2 = df[df['Read'].str.endswith('/2')].add_suffix('_2')
    merged_df = pd.merge(df_1, df_2, left_on='Main_Read_1', right_on='Main_Read_2', how='outer')
    return merged_df.rename(columns={'Main_Read_1': 'Main_Read'})

def determine_final_lca(row):
    all_bids = [row[f'BID{i}_1'] for i in range(1, 4) if pd.notnull(row[f'BID{i}_1'])] + \
               [row[f'BID{i}_2'] for i in range(1, 4) if pd.notnull(row[f'BID{i}_2'])]
    all_kids = [row[f'KID{i}_1'] for i in range(1, 4) if pd.notnull(row[f'KID{i}_1'])] + \
               [row[f'KID{i}_2'] for i in range(1, 4) if pd.notnull(row[f'KID{i}_2'])]
    all_taxids = all_bids + all_kids

    lcas = []
    for taxid1, taxid2 in itertools.combinations(all_taxids, 2):
        if taxid1 and taxid2:
            try:
                lca = ncbi.get_topology([taxid1, taxid2])
                lcas.append((taxid1, taxid2, lca.taxid))
            except KeyError:
                pass

    species_lcas = [(lca, get_rank(lca[2])) for lca in lcas if get_rank(lca[2]) == 'species']
    species_lcas.sort(key=lambda x: (x[1], not is_classified(x[0][2])))

    if species_lcas:
        return species_lcas[0][0]
    elif lcas:
        lcas.sort(key=lambda x: ncbi.get_rank([x[2]])[x[2]])
        return lcas[0]
    else:
        return None, None, None

def fetch_lineage(taxid):
    try:
        lineage = ncbi.get_lineage(taxid)
        return ncbi.get_rank(lineage), ncbi.get_taxid_translator(lineage)
    except ValueError:
        return {}, {}

def build_hierarchy(lca_counts):
    parent_of = {}
    taxon_count = defaultdict(int)
    for lca, count in lca_counts.items():
        try:
            lineage = ncbi.get_lineage(lca)
            prev_taxid = None
            for taxid in lineage:
                if prev_taxid is not None:
                    parent_of[taxid] = prev_taxid
                prev_taxid = taxid
                taxon_count[taxid] += count
        except ValueError:
            print(f'Error processing LCA: {lca}')
    return parent_of, taxon_count

def print_hierarchy_to_file(parent_of, taxon_count, taxid, file, indent=0):
    children = [child for child, parent in parent_of.items() if parent == taxid]
    own_count = taxon_count[taxid] - sum(taxon_count[child] for child in children)

    if own_count > 0 and taxid != 1:
        taxon_name = ncbi.get_taxid_translator([taxid]).get(taxid, 'Unknown')
        file.write(' ' * indent + f'{taxon_name} (TaxID: {taxid}) {own_count}\n')

    for child in children:
        print_hierarchy_to_file(parent_of, taxon_count, child, file, indent + 4)

def translate_taxids_to_names(taxids):
    valid_taxids = [int(taxid) for taxid in taxids if pd.notnull(taxid)]
    names_dict = ncbi.get_taxid_translator(valid_taxids)
    return [names_dict.get(int(taxid), 'Unknown') if pd.notnull(taxid)
    return [names_dict.get(int(taxid), 'Unknown') if pd.notnull(taxid) else 'Unknown' for taxid in taxids]

def main(dataframe, prefix):
    preprocessed_df = preprocess_dataframe(dataframe)
    with ProcessPoolExecutor(max_workers=8) as executor:
        results = list(executor.map(process_row, [row for _, row in preprocessed_df.iterrows()]))

    results_df = pd.DataFrame(results, columns=['Main_Read', 'BID_Used', 'KID_Used', 'LCA'])
    results_df.to_csv(f'{output_dir}/{prefix}_lca_numbers.csv', sep=',', index=False)

    lca_counts = results_df['LCA'].value_counts().to_dict()
    parent_of, taxon_count = build_hierarchy(lca_counts)

    with open(f'{output_dir}/{prefix}_taxonomy_hierarchy.txt', 'w') as file:
        roots = set(parent_of.values()) - set(parent_of.keys())
        for root in roots:
            root_name = ncbi.get_taxid_translator([root]).get(root, 'Unknown')
            file.write(f'Root: {root_name} (TaxID: {root})\n')
            print_hierarchy_to_file(parent_of, taxon_count, root, file)

    results_df['BID_Used'] = translate_taxids_to_names(results_df['BID_Used'])
    results_df['KID_Used'] = translate_taxids_to_names(results_df['KID_Used'])
    results_df['LCA'] = translate_taxids_to_names(results_df['LCA'])

    results_df.to_csv(f'{output_dir}/{prefix}_lca_results.csv', sep=',', index=False)

    return results_df

for _, row in manifest.iterrows():
    prefix = row['prefix']
    print(f'Processing {prefix}')
    dataframe = pd.read_csv(f'process/{prefix}/{prefix}_krakBlastConfReads.csv', sep=',', quotechar='"')
    final_combined_df = main(dataframe, prefix)

print('Taxonomic analysis complete.')
"

# 10. Taxonomic analysis and metrics calculation
echo "Step 10: Calculating taxonomic analysis metrics"
run_python_script "
import sys
import pandas as pd
from ete3 import NCBITaxa

manifest_file = '$manifest_file'
input_dir = '$output_dir'
output_dir = '$output_dir'

# Read manifest file
manifest = pd.read_csv(manifest_file, sep='\t', header=None, names=['read1', 'read2', 'prefix'])

# Initialize NCBITaxa
ncbi = NCBITaxa()

def is_match_at_level(taxid, label_taxid, level):
    if pd.isna(taxid) or pd.isna(label_taxid):
        return False
    if taxid == label_taxid:
        return True
    lineage = ncbi.get_lineage(taxid)
    label_lineage = ncbi.get_lineage(label_taxid)
    label_rank = ncbi.get_rank([label_taxid]).get(label_taxid, None)
    if label_rank != level:
        label_lineage = [ancestor for ancestor in label_lineage if ncbi.get_rank([ancestor])[ancestor] == level]
        label_taxid = label_lineage[0] if label_lineage else None
    return label_taxid in lineage

def calculate_precision(df, label_col, pred_col, level):
    matches = df.apply(lambda x: is_match_at_level(x[pred_col], x[label_col], level) if pd.notna(x[pred_col]) else False, axis=1)
    TP = sum(matches)
    FN = sum(~matches & ~pd.isna(df[label_col]))
    return TP / (TP + FN) if TP + FN > 0 else 0

levels = ['species', 'genus', 'family']

for _, row in manifest.iterrows():
    prefix = row['prefix']
    print(f'Processing {prefix}')
    
    # Load data
    krak_blast_conf_reads_df = pd.read_csv(f'{input_dir}/{prefix}_krakBlastConfReads.csv', sep=',', quotechar='"')
    lca_numbers_df = pd.read_csv(f'{input_dir}/{prefix}_lca_numbers.csv')
    read_label_library_df = pd.read_csv(f'{input_dir}/{prefix}_read_label_library.csv')
    
    # Merge DataFrames
    lca_numbers_df['Read'] = lca_numbers_df['Main_Read'].astype(str).apply(lambda x: x.split('/')[0])
    lca_numbers_df_1 = lca_numbers_df.copy()
    lca_numbers_df_2 = lca_numbers_df.copy()
    lca_numbers_df_1['Read'] = lca_numbers_df_1['Read'] + '/1'
    lca_numbers_df_2['Read'] = lca_numbers_df_2['Read'] + '/2'
    lca_numbers_expanded = pd.concat([lca_numbers_df_1, lca_numbers_df_2])
    
    merged_df = read_label_library_df.merge(krak_blast_conf_reads_df, left_on='ReadID', right_on='Read', how='left')
    merged_df = merged_df.merge(lca_numbers_expanded, on='Read', how='left')
    
    final_df = merged_df[['Read', 'Label', 'BID1', 'KID1', 'LCA']].copy()
    final_df.rename(columns={'BID1': 'BLAST', 'KID1': 'Kraken'}, inplace=True)
    
    # Calculate precision for each level
    precision_data = []
    for level in levels:
        kraken_precision = calculate_precision(final_df, 'Label', 'Kraken', level)
        blast_precision = calculate_precision(final_df, 'Label', 'BLAST', level)
        pipeline_precision = calculate_precision(final_df, 'Label', 'LCA', level)
        precision_data.append([level.capitalize(), kraken_precision, blast_precision, pipeline_precision])
    
    # Create a DataFrame for the precision data
    metrics_df = pd.DataFrame(precision_data, columns=['Precision Level', 'Kraken', 'BLAST', 'Pipeline'])
    
    # Write the DataFrame to a CSV file
    metrics_df.to_csv(f'{output_dir}/{prefix}_pipeline_metrics.csv', index=False)

print('Taxonomic analysis and metrics calculation complete.')
"

echo "Pipeline execution completed."
