#' Count Spikes
#'
#' This function takes a fastq file and counts ocurrences of spike sequences from spike file. Generates three types of output:
#' \enumerate{
#' \item Quality control summary
#' \item A list of fastq read ids that contain spikes
#' \item A table containing spikes and counts for each spike.
#' }
#' @param input.fastq Single-end or merged Paired-end fastq file.
#' @param spike.file Tab-separated file containing spike sequences.
#' @param spike.length Any integer. Generally 9 for spike removal and 25 for normalization.
#' @param output.dir Directory in which output files will be written

#   load depdencies
#.libPaths("/mnt/lustre1/CompBio/lib/R/library")
suppressMessages(source("https://bioconductor.org/biocLite.R", echo=FALSE, verbose=FALSE))
suppressMessages(library(ShortRead));
#library(stringr);
#library(Biostrings);
#biocLite("Biostrings")

#count.spikes <- function(input.fastq, spike.file, spike.length, output.count, output.remove){

input.fastq <- commandArgs(trailingOnly=T)[1]
spike.file <- commandArgs(trailingOnly=T)[2]
spike.length <- commandArgs(trailingOnly=T)[3]
output.count <- commandArgs(trailingOnly=T)[4]
output.remove <- commandArgs(trailingOnly=T)[5]

  #   Read in fastq file
  fastq.reads <- readFastq(input.fastq);
  num.fastqs <- length(fastq.reads);
  cat(num.fastqs, " fastq reads to process\n", sep="");

  #   open spike file
  spike.table <- read.csv(spike.file,
                          stringsAsFactors=FALSE,
                          sep=" ");
  #   extract spike sequences into a vector.
  spikes <- spike.table$SPIKE;

  #   Verify that all spikes are same length
  lengths <- lapply(spikes, nchar)
  len.spike <- nchar(spikes[1])
  comparison <- len.spike == lengths
  comparison <- which(comparison == FALSE)
  if (length(comparison) > 0){
    stop("At least one spike is incorrect length")
  }

#   TODO:  add some error-checking, verifying that all spikes are of same length, etc.

  #    trim to length (regardless of whether we're using the spike or its reverse-compl.
  spikes <- strtrim(spikes, spike.length);
  
  #   create outputs
  output.table <- spike.table;
  output.table$spike.count <- 0;
  records.to.remove.ids <- character();

  #   count spikes - note there are two cases
  for(i in 1:length(spikes))  {
     
      #   get a single spike 
      current.spike <- spikes[i];
      #   count spike occurences
      current.spike.counts <- vcountPattern(current.spike, sread(fastq.reads), max.mismatch = 1, with.indels = TRUE);
      #   append to output table
      output.table[i,]$spike.count <- sum(as.logical(current.spike.counts));
      #   update progress
      cat("\t", output.table[i,]$spike.count, " spikes detected\n", sep="");
      #   record readids (as they'll need removed later; we use the readid to do so)
      records.to.remove <- which(current.spike.counts > 0);
      if(length(records.to.remove) > 0)   {
          temp.reads <- fastq.reads[records.to.remove];
          records.to.remove.ids <- c(records.to.remove.ids, 
                                      as.character(temp.reads@id));
      }   #   fi
          
      #   In the case where spike.length == 9, all spikes are identical (since the
      #       first nine characters of each spike are identical.  We use this break
      #       to short-circuit the "for" loop, since we get all the information we
      #       need in the first pass
      #   If the first 9 characters are not identical, stop.
      if(spike.length == 9){
        check <- spikes[1] == spikes
        check <- which(check == FALSE)
        if (length(check) == 0){
          break;
        } else {
          stop("First nine characters of spikes are not identical.")
        }   #   else
      }   #   fi
            
      #   update progress to the user 
      if((i %% 10) == 0)  {
          cat("Processing spike ", i, " out of ", length(spikes), "\n", sep="");
      }   #   fi
  }   #   for i
  
    
  #   build names of output files
  #   strip ".fastq" from file name.  Not necessary, just more hygienic 
#  filename.no.extension <- sub("[.][^.]*$", "", basename(input.fastq));
  #   construct spike count names
#  output.count <- paste(filename.no.extension, ".spike.counts.", spike.length, "bp.txt", sep="");
  #   construct to remove names
#  output.remove <- paste(filename.no.extension, ".reads.to.remove.txt", sep="");

  #   write outputs
  write.table(output.table,
              file=output.count,
              quote=FALSE,
              sep=",");
  
  write.table(records.to.remove.ids,
              file=output.remove,
              quote=FALSE,
              row.names=FALSE);
#              col.names=FALSE);
  

