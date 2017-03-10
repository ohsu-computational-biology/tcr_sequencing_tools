# This script is designed to extract the read IDs that were assigned to the V22 pseudogene
# during the assembly process

# First we take the export_clones file and extract all of the Clone IDs for clones whose
# best V hit is V22.

# Then we take the export_alignment file and extract all of the read IDs for the previously
# extracted clone IDs

# Inputs: export clones file and export alignment file
# Outputs: file listing all read IDs

###
### Load dependencies ###
###

.libPaths("/mnt/lustre1/CompBio/lib/R/library")
library(ShortRead)
library(stringr)

###
### Set command line arguments ###
###

arguments <- commandArgs(trailingOnly = TRUE)

align.file <- arguments[1]
assemble.file <- arguments[2]
fastq.file <- arguments[3]
output.dir <- arguments[4]

###
### Read in data ###
###

align <- read.delim(align.file, sep = '\t', header = T, na.strings = c('', ' '), 
                    stringsAsFactors = F)


assemble <- read.delim(assemble.file, sep = '\t', header = T, na.strings = c('', ' '), 
                       stringsAsFactors = F)

fastq <- readFastq(fastq.file)

###
### Metadata and Error Checking ###
###

# Extract Sample Numbers
align.sample.number <- gsub(".*_S|_alignment.*", '', align.file)

assemble.sample.number <- gsub(".*_S|_alignment.*", '', assemble.file)

fastq.sample.number <- gsub(".*_S|\\..*", '', fastq.file)

# Error Check
if (align.sample.number != assemble.sample.number) {
   stop("Alignment and Assemble files do not match.")
} # if

if (align.sample.number != fastq.sample.number) {
   stop("Alignment and Fastq files do not match.")
} # if

if (assemble.sample.number != fastq.sample.number) {
   stop("Assemble and Fastq files do not match.")
} # if

sample.number <- align.sample.number


###
### Pre-process ###
###

# Clean up assemble columns
assemble$Best.V.hit <- sub("TRB", "", assemble$Best.V.hit)
assemble$Best.V.hit <- sub("\\*00", "", assemble$Best.V.hit)

###
### Obtain V22 Read IDs ###
###

# Extract clone ids that correspond to V22
clone.ids <- assemble[(assemble$Best.V.hit == "V22"), 1]

# Extract Read ids from corresponding clone ids
read.ids <- align[(align$Clone.Id %in% clone.ids), "Description.R1"]

###
### Extract matching reads from fastq file 
###

if (length(read.ids) > 0){
   id.capture <- srFilter(function(x){
   	     (x@id %in% read.ids)
	     }, # function(x)
	     name = "id.capture");
   output.fastq <- fastq[id.capture(fastq)]
} else {
   output.fastq <- NULL
   print("There were no clones assigned to V22")
} # if

###
### Create output file name and write to file ###
###

output.name <- paste("S", sample.number, "_V22_reads.fastq", sep = '')

if (is.null(output.fastq)) {
   cat("No reads assigned to V22, so no output file", '\n')
} else {
  cat("Writing output to: ", output.dir, output.name, '\n', sep = '')

  writeFastq(output.fastq, paste(output.dir, output.name, sep = ''), compress = FALSE)
} # if
