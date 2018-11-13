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
keep_v <- arguments[4]     # logical. TRUE - keep the reads specified in read_dir_v files. FALSE - remove those reads and keep the rest.

### Default to removing reads
if (is.na(keep_v)) keep_v <- FALSE

### Get files
fastq_files_v <- list.files(fastq_dir_v)
read_files_v <- list.files(read_dir_v)

### subset fastq files to only contain reads samples
### Normally doesn't do anything, but just in case there are different numbers of samples
samples_v <- sapply(read_files_v, function(x) gsub("^.*_S|_read_ids.txt|_contamIDs.txt", "", x), USE.NAMES = F)
fastq_files_v <- fastq_files_v[gsub("^.*_S|\\.assemb.*|_R1.*|_R2.*", "", fastq_files_v) %in% samples_v]

### Sort by sample number
fastq_files_v <- fastq_files_v[order(as.numeric(gsub("^.*_S|\\.assemb.*|_R1.*|_R2.*", "", fastq_files_v)))]
read_files_v <- read_files_v[order(as.numeric(gsub("^.*_S|_read_ids.txt|_contamIDs.txt", "", read_files_v)))]

### Split forward and reverse reads
fwd_files_v <- grep("_R1_", fastq_files_v, value = T)
rev_files_v <- grep("_R2_", fastq_files_v, value = T)

### Extract reads from each fastq file using corresponding read ID file
for (i in 1:length(read_files_v)){
  ## Assemble output name
  name_v <- unlist(strsplit(read_files_v[i], split = "_"))
  fwd_out_name_v <- paste(name_v[1], name_v[2], "R1_decontam.fastq", sep = "_")
  rev_out_name_v <- paste(name_v[1], name_v[2], "R2_decontam.fastq", sep = "_")
  
  ## Get fastq data and assemble read IDs into a character vector
  fwd.fastq.records <- readFastq(paste0(fastq_dir_v, fwd_files_v[i]))
  rev.fastq.records <- readFastq(paste0(fastq_dir_v, rev_files_v[i]))
  ids.to.remove <- unlist(fread(paste0(read_dir_v, read_files_v[i]), sep = '\t')[,1, with = F], use.names = F)

  alt.ids.to.remove <- gsub(" 1", " 2", ids.to.remove)
  ids.to.remove <- c(ids.to.remove, alt.ids.to.remove)
  
  ## Define function to filter out reads for export
  if (keep_v) {
    id.filter <- srFilter(function(x) {
      (x@id %in% ids.to.remove)
    },
    name="id.filter")
  } else {
    id.filter <- srFilter(function(x) {
      (!(x@id %in% ids.to.remove))
    },
    name="id.filter")
  } # fi

  ## Apply filtering function
  fwd.output.fastq.records <- fwd.fastq.records[id.filter(fwd.fastq.records)]
  rev.output.fastq.records <- rev.fastq.records[id.filter(rev.fastq.records)]
  
  ## Write output
  writeFastq(fwd.output.fastq.records, paste0(out_dir_v, fwd_out_name_v), compress = F)
  writeFastq(rev.output.fastq.records, paste0(out_dir_v, rev_out_name_v), compress = F)
} # for

