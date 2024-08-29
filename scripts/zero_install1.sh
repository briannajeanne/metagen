#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Update package lists
sudo apt-get update

# Install basic dependencies
sudo apt-get install -y \
    python3 \
    python3-pip \
    r-base \
    gzip

# Install InSilicoSeq
pip3 install InSilicoSeq

# Install required R packages
sudo Rscript -e 'install.packages(c("optparse", "readr", "dplyr"), repos="https://cloud.r-project.org")'

# Verify installations
echo "Verifying installations..."

python3 --version
pip3 --version
R --version
iss --version
gzip --version

echo "Installation complete. Please ensure that all commands above displayed version information."
echo "If any command was not found, please check the installation logs and try again."
