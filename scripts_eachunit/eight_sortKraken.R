#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript sortKraken.R <manifest_file> <output_dir>")
}

manifest_file <- args[1]
output_dir <- args[2]

manifest <- read.delim(manifest_file, header = FALSE, col.names = c("read1", "read2", "prefix"))

for (i in 1:nrow(manifest)) {
  prefix <- manifest$prefix[i]
  setwd(file.path("process", prefix))
  
  sampleName <- prefix
  
  # Rest of the script remains the same, just replace hardcoded filenames with prefixed versions
  # For example:
  # read.delim('krakUniq_sample.report', comment.char = "#")
  # becomes:
  # read.delim(paste0(prefix, '_krakUniq_sample.report'), comment.char = "#")
  
  # At the end, use prefixed output filenames:
  write.table(allKrakenReads, paste0(prefix, '_krakenSelVirReads.tsv'), row.names = F, sep = '\t', quote = F)
  write.table(reads, paste0(prefix, '_krakenReads.id'), row.names = F, col.names = F, quote = F)
  write_csv(result, paste0(prefix, '_krakenSelVirReads.tsv'))
  
  setwd("../..")
}
