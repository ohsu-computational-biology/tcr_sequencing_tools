#	Wrapper script for optimized.remove.spikes.R
#
#	Loads optimized.remove.spikes.R, then calls remove.fastqs() for all files in directory
#
#	This script should be placed in the same directory containing, for each sample SX, files:
#
#		SX.fastq
#		SX.reads.to.remove.txt
#
#	Load required libraries
library(stringr);

merged.reads <- TRUE;

#	Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);
working.dir <- arguments[1];

#	Specify relevant details
source("/mnt/lustre1/CompBio/genomic_resources/tcrseq/optimized.remove.spikes.R");
fastq.suffix <- ".fastq";
reads.to.remove.suffix <- ".reads.to.remove.txt";

#	Examine the current directory for the files to process
files.in.dir <- list.files(working.dir);

#	Prepare a list of files to process
#   text file indicating reads to remove
reads.to.remove.files <- str_detect(files.in.dir, reads.to.remove.suffix);
reads.to.remove.files <- files.in.dir[reads.to.remove.files];
reads.to.remove.files <- paste(working.dir, reads.to.remove.files, sep="");
#   Untouched fastq files
fastq.files <- str_detect(files.in.dir, fastq.suffix);
fastq.files <- files.in.dir[fastq.files];
fastq.files <- paste(working.dir, fastq.files, sep="");
file.roots.1 <- unique(reads.to.remove.files);
file.roots.2 <- unique(fastq.files);

#	error-check
if(length(file.roots.1) != length(file.roots.2))
	stop("Mismatch between number of fastq files and number of reads.to.remove.files in directory")

for(i in 1:length(fastq.files))	{
    #   get a fastq file to process
    curr.fastq <- fastq.files[i];
  
    #   get the appropriate read.to.remove files 
     if(merged.reads)    {
        curr.reads.forward <- reads.to.remove.files[i];
        curr.reads.reverse <- reads.to.remove.files[i];
    }   else {  #   for each fastq file there are two reads.to.remove.files, one forward, one reverse
        #	list.files() returns files in alphabetical order.  We process three files at a time,
        #		assuming this order (though also confirming the order)
        if((i %% 2) == 0)	{
            curr.reads.forward <- reads.to.remove.files[i - 1];
            curr.reads.reverse <- reads.to.remove.files[i];
        }	else	{
            curr.reads.forward <- reads.to.remove.files[i];
            curr.reads.reverse <- reads.to.remove.files[i + 1];
        }	#	else
    }   #   esle

	#	call function
	remove.fastqs(curr.fastq, curr.reads.forward, curr.reads.reverse);

	if((i %% 10) == 0)	{
		cat("Processed: ", i, " (of ", length(fastq.files), ")\n",
			"\t", curr.fastq, "\n",
			"\t", curr.reads.forward, "\n",
			"\t", curr.reads.reverse, "\n\n",
			sep="");
	}	#	fi

}	#	for i


