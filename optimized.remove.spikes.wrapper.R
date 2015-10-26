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
read.files <- str_detect(files.in.dir, reads.to.remove.suffix);
read.files <- files.in.dir[read.files];
read.files <- paste(working.dir, read.files, sep="");
fastq.files <- str_detect(files.in.dir, fastq.suffix);
fastq.files <- files.in.dir[fastq.files];
fastq.files <- paste(working.dir, fastq.files, sep="");
file.roots.1 <- unique(read.files);
file.roots.2 <- unique(fastq.files);

#	error-check
if(length(file.roots.1) != length(file.roots.2))
	stop("Mismatch between number of fastq files and number of read files in directory")

for(i in 1:length(fastq.files))	{

	#	list.files() returns files in alphabetical order.  We process three files at a time,
	#		assuming this order (though also confirming the order)
	curr.fastq <- fastq.files[i];
	if((i %% 2) == 0)	{
		curr.reads.forward <- read.files[i - 1];
		curr.reads.reverse <- read.files[i];
	}	else	{
		curr.reads.forward <- read.files[i];
		curr.reads.reverse <- read.files[i + 1];
	}	#	else

	#	call function
	remove.fastqs(curr.fastq, curr.reads.forward, curr.reads.reverse);

	if((i %% 10) == 0)	{
		cat("Processed:\n",
			"\t", curr.fastq, "\n",
			"\t", curr.reads.forward, "\n",
			"\t", curr.reads.reverse, "\n\n");
	}	#	fi

}	#	for i


