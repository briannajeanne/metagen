#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Pipeline parameter defaults
params.manifest = "path/to/manifest.tsv"
params.adapter_fasta = "path/to/adapter.fasta"
params.krakenuniq_db = "path/to/krakenuniq_db"
params.blast_db = "path/to/blast_db"
params.output_dir = "results"

// Define input channels
Channel
    .fromPath(params.manifest)
    .splitCsv(header:false, sep:'\t')
    .map{ row -> tuple(row[2], file(row[0]), file(row[1])) }
    .set { input_samples }

// Process 1: Filter and merge reads
process filterAndMerge {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id), path("${sample_id}_merged_reads.fastq"), emit: merged_reads
    
    script:
    """
    fastp -i $read1 -I $read2 \
    --detect_adapter_for_pe \
    --adapter_fasta ${params.adapter_fasta} \
    --merge --merged_out ${sample_id}_merged_reads.fastq \
    --include_unmerged \
    -l 100 \
    --thread ${task.cpus}
    """
}

// Process 2: Convert FASTQ to FASTA
process fastqToFasta {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(merged_fastq)
    
    output:
    tuple val(sample_id), path("${sample_id}_merged_reads.fa"), emit: merged_fasta
    
    script:
    """
    seqtk seq -a $merged_fastq > ${sample_id}_merged_reads.fa
    """
}

// Process 3: Run KrakenUniq
process runKrakenUniq {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(merged_fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}_krakUniq_sample.report"), path("${sample_id}_krakUniq_sample.kraken"), path("${sample_id}_krakUniq_classified.reads"), emit: kraken_output
    
    script:
    """
    krakenuniq --threads ${task.cpus} --db ${params.krakenuniq_db} \
    --report-file ${sample_id}_krakUniq_sample.report --output ${sample_id}_krakUniq_sample.kraken \
    --classified-out ${sample_id}_krakUniq_classified.reads $merged_fasta
    """
}

// Process 4: Sort Kraken results
process sortKraken {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(kraken_report), path(kraken_output), path(classified_reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_krakenSelVirReads.tsv"), path("${sample_id}_krakenReads.id"), emit: sorted_kraken
    
    script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    
    sampleName <- "$sample_id"
    
    virusId <- function() {
        df <- read.delim("$kraken_report", comment.char = "#")
        virus <- as.data.frame(df)
        colnames(virus)[colnames(virus) == "taxName"] <- "Virus"
        virusID <- data.frame(unique(virus[, c('taxID', 'Virus')]))
        colnames(virusID)[1] <- 'krakTax'
        return(virusID)
    }
    
    idList <- virusId()
    readNameTab <- read.delim("$kraken_output", F, sep = '\t')
    colnames(readNameTab)[3] <- 'krakTax'
    colnames(readNameTab)[2] <- 'Read'
    selCol <- readNameTab[, c('krakTax', 'Read')]
    cReads <- plyr::join(idList, selCol, by = 'krakTax', type = 'left', match = 'all')
    allKrakenReads <- unique(cReads[, c('Virus', 'Read')])
    allKrakenReads\$Sample <- sampleName
    
    write.table(allKrakenReads, '${sample_id}_krakenSelVirReads.tsv', row.names = F, sep = '\t', quote = F)
    
    reads <- unique(allKrakenReads[, 'Read'])
    write.table(reads, '${sample_id}_krakenReads.id', row.names = F, col.names = F, quote = F)
    
    file_path <- "$kraken_output"
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
    write_csv(result, '${sample_id}_krakenSelVirReads.tsv')
    """
}

// Process 5: Extract Kraken classified reads
process extractKrakenReads {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(kraken_reads), path(kraken_ids), path(classified_reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_selectKraken.reads"), emit: selected_reads
    
    script:
    """
    seqtk subseq $classified_reads $kraken_ids > ${sample_id}_selectKraken.reads
    """
}

// Process 6: Split reads
process splitReads {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(selected_reads)
    
    output:
    tuple val(sample_id), path("splitSeq10K"), emit: split_reads
    
    script:
    """
    mkdir splitSeq10K
    seqkit split $selected_reads -s 900 -j ${task.cpus} -O ./splitSeq10K/
    cd splitSeq10K
    for j in *.reads; do mkdir "\${j}_dir"; mv "\$j" "./\${j}_dir/reads.fas"; done
    """
}

// Process 7: Run BLAST
process runBlast {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(split_reads_dir)
    
    output:
    tuple val(sample_id), path("${sample_id}_blast_results"), emit: blast_results
    
    script:
    """
    mkdir ${sample_id}_blast_results
    for dir in $split_reads_dir/*/
    do
        blastn -db ${params.blast_db} -query \$dir/reads.fas \
        -num_threads ${task.cpus} -out ${sample_id}_blast_results/\$(basename \$dir)_blast.results \
        -outfmt '7 staxids qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 10
    done
    """
}

// Process 8: Taxonomic analysis
process taxonomicAnalysis {
    tag "$sample_id"
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(blast_results), path(kraken_results)
    
    output:
    tuple val(sample_id), path("${sample_id}_lca_results.csv"), path("${sample_id}_taxonomy_hierarchy.txt"), emit: tax_results
    
    script:
    """
    #!/usr/bin/env python3
    import os
    import pandas as pd
    from ete3 import NCBITaxa
    import itertools
    from concurrent.futures import ProcessPoolExecutor
    from collections import defaultdict

    ncbi = NCBITaxa()

    def get_rank(taxid):
        return ncbi.get_rank([taxid]).get(taxid, None)

    def is_classified(taxid):
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        return "unclassified" not in " ".join(names.values()).lower()

    def get_lca(taxids):
        try:
            combinations = list(itertools.combinations(taxids, 2))
            lcas = [(pair, ncbi.get_topology(pair)) for pair in combinations if pair[0] and pair[1]]
            if not lcas:
                return None, None, None
            lcas.sort(key=lambda x: ncbi.get_rank([x[1].taxid])[x[1].taxid])
            return lcas[0][0][0], lcas[0][0][1], lcas[0][1].taxid
        except KeyError as e:
            print(f"KeyError encountered: {e}")
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
        all_bids = [row[f"BID{i}_1"] for i in range(1, 4) if pd.notnull(row[f"BID{i}_1"])] + \
                   [row[f"BID{i}_2"] for i in range(1, 4) if pd.notnull(row[f"BID{i}_2"])]
        all_kids = [row[f"KID{i}_1"] for i in range(1, 4) if pd.notnull(row[f"KID{i}_1"])] + \
                   [row[f"KID{i}_2"] for i in range(1, 4) if pd.notnull(row[f"KID{i}_2"])]
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
                print(f"Error processing LCA: {lca}")
        return parent_of, taxon_count

    def print_hierarchy_to_file(parent_of, taxon_count, taxid, file, indent=0):
        children = [child for child, parent in parent_of.items() if parent == taxid]
        own_count = taxon_count[taxid] - sum(taxon_count[child] for child in children)

        if own_count > 0 and taxid != 1:
            taxon_name = ncbi.get_taxid_translator([taxid]).get(taxid, "Unknown")
            file.write(" " * indent + f"{taxon_name} (TaxID: {taxid}) {own_count}\\n")

        for child in children:
            print_hierarchy_to_file(parent_of, taxon_count, child, file, indent + 4)

    def translate_taxids_to_names(taxids):
        valid_taxids = [int(taxid) for taxid in taxids if pd.notnull(taxid)]
        names_dict = ncbi.get_taxid_translator(valid_taxids)
        return [names_dict.get(int(taxid), "Unknown") if pd.notnull(taxid) else "Unknown" for taxid in taxids]

    def main(dataframe, prefix):
        preprocessed_df = preprocess_dataframe(dataframe)
        with ProcessPoolExecutor(max_workers=8)
with ProcessPoolExecutor(max_workers=8) as executor:
            results = list(executor.map(process_row, [row for _, row in preprocessed_df.iterrows()]))

        results_df = pd.DataFrame(results, columns=['Main_Read', 'BID_Used', 'KID_Used', 'LCA'])
        results_df.to_csv(f'{prefix}_lca_numbers.csv', sep=',', index=False)

        lca_counts = results_df['LCA'].value_counts().to_dict()
        parent_of, taxon_count = build_hierarchy(lca_counts)

        with open(f'{prefix}_taxonomy_hierarchy.txt', 'w') as file:
            roots = set(parent_of.values()) - set(parent_of.keys())
            for root in roots:
                root_name = ncbi.get_taxid_translator([root]).get(root, "Unknown")
                file.write(f"Root: {root_name} (TaxID: {root})\\n")
                print_hierarchy_to_file(parent_of, taxon_count, root, file)

        results_df['BID_Used'] = translate_taxids_to_names(results_df['BID_Used'])
        results_df['KID_Used'] = translate_taxids_to_names(results_df['KID_Used'])
        results_df['LCA'] = translate_taxids_to_names(results_df['LCA'])

        results_df.to_csv(f'{prefix}_lca_results.csv', sep=',', index=False)

        return results_df

    if __name__ == "__main__":
        dataframe = pd.read_csv("${sample_id}_krakBlastConfReads.csv", sep=',', quotechar='"')
        final_combined_df = main(dataframe, "${sample_id}")
        print("Taxonomic analysis complete.")
    """
}

// Process 9: Calculate taxonomic metrics
process calculateMetrics {
    tag "$sample_id"
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(lca_results), path(taxonomy_hierarchy)
    
    output:
    path "${sample_id}_pipeline_metrics.csv"
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    from ete3 import NCBITaxa

    def is_match_at_level(taxid, label_taxid, level):
        if pd.isna(taxid) or pd.isna(label_taxid):
            return False
        if taxid == label_taxid:
            return True
        ncbi = NCBITaxa()
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

    def calculate_metrics(sample_id):
        # Load data
        krak_blast_conf_reads_df = pd.read_csv(f'{sample_id}_krakBlastConfReads.csv', sep=',', quotechar='"')
        lca_numbers_df = pd.read_csv(f'{sample_id}_lca_numbers.csv')
        read_label_library_df = pd.read_csv(f'{sample_id}_read_label_library.csv')
        
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
        levels = ['species', 'genus', 'family']
        precision_data = []
        for level in levels:
            kraken_precision = calculate_precision(final_df, 'Label', 'Kraken', level)
            blast_precision = calculate_precision(final_df, 'Label', 'BLAST', level)
            pipeline_precision = calculate_precision(final_df, 'Label', 'LCA', level)
            precision_data.append([level.capitalize(), kraken_precision, blast_precision, pipeline_precision])
        
        # Create a DataFrame for the precision data
        metrics_df = pd.DataFrame(precision_data, columns=['Precision Level', 'Kraken', 'BLAST', 'Pipeline'])
        
        # Write the DataFrame to a CSV file
        metrics_df.to_csv(f'{sample_id}_pipeline_metrics.csv', index=False)

    if __name__ == "__main__":
        calculate_metrics("${sample_id}")
        print("Metrics calculation complete.")
    """
}

// Main workflow
workflow {
    filterAndMerge(input_samples)
    fastqToFasta(filterAndMerge.out.merged_reads)
    runKrakenUniq(fastqToFasta.out.merged_fasta)
    sortKraken(runKrakenUniq.out.kraken_output)
    extractKrakenReads(sortKraken.out.sorted_kraken.join(runKrakenUniq.out.kraken_output))
    splitReads(extractKrakenReads.out.selected_reads)
    runBlast(splitReads.out.split_reads)
    taxonomicAnalysis(runBlast.out.blast_results.join(sortKraken.out.sorted_kraken))
    calculateMetrics(taxonomicAnalysis.out.tax_results)
}
