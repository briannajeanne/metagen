#!/bin/bash

set -e

# Define installation directory
INSTALL_DIR="$HOME/bioinfo_pipeline"
mkdir -p $INSTALL_DIR
cd $INSTALL_DIR

# Install Miniconda if not already installed
if ! command -v conda &> /dev/null; then
    echo "Installing Miniconda..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $HOME/miniconda
    echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> $HOME/.bashrc
    source $HOME/.bashrc
fi

# Create conda environment
conda create -n bioinfo_pipeline python=3.8 -y
conda activate bioinfo_pipeline

# Install software using conda
conda install -y -c bioconda \
    fastp \
    seqtk \
    blast \
    seqkit \
    ete3 \
    pandas \
    r-base \
    r-tidyverse \
    r-plyr

# Install KrakenUniq
echo "Installing KrakenUniq..."
conda install -y -c bioconda krakenuniq

# Install Nextflow
echo "Installing Nextflow..."
conda install -y -c bioconda nextflow

# Create directories for databases
mkdir -p $INSTALL_DIR/databases

# Download and build Kraken database
echo "Downloading and building Kraken database..."
mkdir -p $INSTALL_DIR/databases/kraken_db
cd $INSTALL_DIR/databases/kraken_db
krakenuniq-download --db $INSTALL_DIR/databases/kraken_db standard
krakenuniq-build --db $INSTALL_DIR/databases/kraken_db --threads 4

# Download BLAST nt database
echo "Downloading BLAST nt database..."
mkdir -p $INSTALL_DIR/databases/blast_db
cd $INSTALL_DIR/databases/blast_db
update_blastdb.pl --decompress nt

# Download adapter sequences
echo "Downloading adapter sequences..."
mkdir -p $INSTALL_DIR/databases/adapters
cd $INSTALL_DIR/databases/adapters
wget https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/NexteraPE-PE.fa

# Create a configuration file
echo "Creating configuration file..."
cat > $INSTALL_DIR/pipeline_config.sh << EOF
export BIOINFO_PIPELINE_DIR="$INSTALL_DIR"
export KRAKEN_DB="$INSTALL_DIR/databases/kraken_db"
export BLAST_DB="$INSTALL_DIR/databases/blast_db/nt"
export ADAPTER_FASTA="$INSTALL_DIR/databases/adapters/NexteraPE-PE.fa"
EOF

echo "Installation complete. Please add the following line to your .bashrc file:"
echo "source $INSTALL_DIR/pipeline_config.sh"

echo "Then, run 'source ~/.bashrc' to update your current session."

echo "To run the pipeline, activate the conda environment with:"
echo "conda activate bioinfo_pipeline"

echo "And run Nextflow with:"
echo "nextflow run main.nf --manifest path/to/your/manifest.tsv --adapter_fasta \$ADAPTER_FASTA --krakenuniq_db \$KRAKEN_DB --blast_db \$BLAST_DB --output_dir results"
