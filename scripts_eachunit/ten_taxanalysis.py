import sys
import pandas as pd
from ete3 import NCBITaxa

if len(sys.argv) != 4:
    print("Usage: python taxanalysis.py <manifest_file> <input_dir> <output_dir>")
    sys.exit(1)

manifest_file = sys.argv[1]
input_dir = sys.argv[2]
output_dir = sys.argv[3]

# Read manifest file
manifest = pd.read_csv(manifest_file, sep='\t', header=None, names=['read1', 'read2', 'prefix'])

# Initialize NCBITaxa
ncbi = NCBITaxa()

# Rest of the script remains largely the same, but we'll process each sample separately:

for _, row in manifest.iterrows():
    prefix = row['prefix']
    
    # Load data
    krak_blast_conf_reads_df = pd.read_csv(f"{input_dir}/{prefix}_krakBlastConfReads.csv", sep=',', quotechar='"')
    lca_numbers_df = pd.read_csv(f"{input_dir}/{prefix}_lca_numbers.csv")
    read_label_library_df = pd.read_csv(f"{input_dir}/{prefix}_read_label_library.csv")
    
    # Process the data (main function logic here)
    
    # Write output files with prefixes
    metrics_df.to_csv(f"{output_dir}/{prefix}_pipeline_metrics.csv", index=False)

print("Taxonomic analysis and metrics calculation complete.")
