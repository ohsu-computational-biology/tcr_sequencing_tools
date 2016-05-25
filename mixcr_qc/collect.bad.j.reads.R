# This script will extract read IDs that aligned in MiXCR from the fastq file so that
# a new fastq with only failed reads is produced.

# Inputs:
#  1. S##_read_ids.txt
#        These are produced by the mixcr.qc.R script
#        This script uses the aligned reads, not the assembled
#  2. Corresponding fastq file

# Outputs:
#  1. fastq file with reads removed
#  2. Summary of this script

# load dependencies
.libPaths("/mnt/lustre1/CompBio/lib/R/library")
library(ShortRead);
library(stringr);

# Set commandline arguments

arguments <- commandArgs(trailingOnly = TRUE);

input.fastq <- arguments[1]
reads.to.remove <- arguments[2]
output.dir <- arguments[3]

# For testing
input.fastq <- "~/Desktop/OHSU/tcr_spike/data/equiv_DNA160107LC/remove_test/DNA160107LC_S1.assembled.fastq"
reads.to.remove <- "~/Desktop/OHSU/tcr_spike/data/equiv_DNA160107LC/remove_test/S1_align_read_ids.txt"
output.dir <- "~/Desktop/OHSU/tcr_spike/data/equiv_DNA160107LC/remove_test/"

# Read in fastq file
fastq.reads <- readFastq(input.fastq)
cat(length(fastq.reads), " fastq reads to process\n", sep = '')

# Create variable for IDs to remove
ids.to.remove <- NULL

# Check file size of reads to remove, and if so, generate alternative IDs
file.size.to.remove <- file.size(reads.to.remove)
if (file.size.to.remove > 0){
  ids.to.remove <- read.delim(reads.to.remove, header = F, stringsAsFactors = F)
  }

# Remove from fastq file
if (length(ids.to.remove$V1) > 0){
    id.filter <- srFilter(function(x){
      !(x@id %in% ids.to.remove$V1)
    }, # function(x)
    name = "id.filter");
  output.fastq.reads <- fastq.reads[id.filter(fastq.reads)]
} else {
  output.fastq.reads <- fastq.reads
} # else

# Output
cat("Extracting ", length(output.fastq.reads), " reads with failed J alignments\n",
    sep = '')

base.name <- strsplit(input.fastq, split = '/')[[1]]
base.name <- strsplit(base.name[length(base.name)], split = '\\.')[[1]][1]
output.name <- paste(base.name, ".failed.j.fastq", sep = '')

cat("Writing output to: ", output.dir, output.name, '\n', sep = '')

writeFastq(output.fastq.reads, paste(output.dir, output.name, sep = ''),
           compress=FALSE)