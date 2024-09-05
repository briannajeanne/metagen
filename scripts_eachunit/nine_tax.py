import os
import sys
from ete3 import NCBITaxa
import pandas as pd
import itertools
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict

if len(sys.argv) != 3:
    print("Usage: python tax.py <manifest_file> <output_dir>")
    sys.exit(1)

manifest_file = sys.argv[1]
output_dir = sys.argv[2]

# Read manifest file
manifest = pd.read_csv(manifest_file, sep='\t', header=None, names=['read1', 'read2', 'prefix'])

# Initialize NCBITaxa
ncbi = NCBITaxa()

# Rest of the script remains largely the same, but we'll process each sample separately:

for _, row in manifest.iterrows():
    prefix = row['prefix']
    
    # Read input file
    input_file = f"process/{prefix}/{prefix}_krakBlastConfReads.csv"
    dataframe = pd.read_csv(input_file, sep=',', quotechar='"')
    
    # Process the dataframe (main function logic here)
    
    # Write output files with prefixes
    results_df.to_csv(f"{output_dir}/{prefix}_lca_numbers.csv", sep=',', index=False)
    
    with open(f"{output_dir}/{prefix}_taxonomy_hierarchy.txt", 'w') as file:
        # Write hierarchy
    
    results_df.to_csv(f"{output_dir}/{prefix}_lca_results.csv", sep=',', index=False)

print("Taxonomic analysis complete.")
