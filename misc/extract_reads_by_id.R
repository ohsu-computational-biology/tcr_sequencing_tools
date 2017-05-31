###
### Extract reads by IDs
###

### Given a fastq file and a one-column text file of read IDs, extract from the fastq file all of the corresponding read information for the reads
### in the read ID file. Output into a fastq file containing only those reads specified in the read IDs file.

### load depdencies
library(ShortRead);
library(stringr);
library(data.table)

### Suppress warnings
options(warn=-1);


### Arguments
arguments <- commandArgs(trailingOnly = T)
fastq_dir_v <- arguments[1]
read_dir_v <- arguments[2]
out_dir_v <- arguments[3]
#fastq_dir_v <- "~/Desktop/OHSU/tcr_spike/data/LIB170111LC/despiked/"
#read_dir_v <- "~/Desktop/OHSU/tcr_spike/data/LIB170111LC/eval/"
#out_dir_v <- "~/Desktop/OHSU/tcr_spike/data/LIB170111LC/assembled_reads/"

### Get files
fastq_files_v <- list.files(fastq_dir_v)
read_files_v <- list.files(read_dir_v)

### subset fastq files to only contain reads samples
samples_v <- sapply(read_files_v, function(x) gsub("^.*_S|_read_ids.txt", "", x), USE.NAMES = F)
fastq_files_v <- fastq_files_v[gsub("^.*_S|\\.assemb.*", "", fastq_files_v) %in% samples_v]

### Sort by sample number
fastq_files_v <- fastq_files_v[order(as.numeric(gsub("^.*_S|\\.assemb.*", "", fastq_files_v)))]
read_files_v <- read_files_v[order(as.numeric(gsub("^.*_S|_read_ids.txt", "", read_files_v)))]


### Extract reads from each fastq file using corresponding read ID file
for (i in 1:length(fastq_files_v)){
  ## Assemble output name
  name_v <- unlist(strsplit(read_files_v[i], split = "_"))
  out_name_v <- paste(name_v[1], name_v[2], "assembled_reads.fastq", sep = "_")
  
  ## Get fastq data and assemble read IDs into a character vector
  fastq.records <- readFastq(paste0(fastq_dir_v, fastq_files_v[i]))
  ids.to.remove <- unlist(fread(paste0(read_dir_v, read_files_v[i]), sep = '\t')[,1, with = F], use.names = F)
  alt.ids.to.remove <- gsub(" 1", " 2", ids.to.remove)
  ids.to.remove <- c(ids.to.remove, alt.ids.to.remove)
  
  ## Define function to filter out reads for export
  id.filter <- srFilter(function(x) {
    (x@id %in% ids.to.remove)
  },
  name="id.filter")
  
  ## Apply filtering function
  output.fastq.records <- fastq.records[id.filter(fastq.records)]
  
  ## Write output
  writeFastq(output.fastq.records, paste0(out_dir_v, out_name_v), compress = F)
} # for

