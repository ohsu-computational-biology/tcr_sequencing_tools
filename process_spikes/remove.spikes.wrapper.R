#	Wrapper script for optimized.remove.spikes.R
#
#	Loads optimized.remove.spikes.R, then calls remove.fastqs() (which is a function defined
#       in that file) for all files in directory specified as an input to this Rscript
#
#	This script should be placed in the same directory containing, for each sample SX, files:
#
#		SX.fastq
#		SX.reads.to.remove.txt
#
#   This script expects as input the path to this directory
#   
#   There are two modes of operation, depending on whether or not merged.reads is TRUE or FALSE.
#   The original workflow kept forward and reverse reads separate.  For example, for a single
#       sample:
#
#           S1_R1.fastq
#           S1_R2.fastq
#           S1_R1.reads.to.remove.txt
#           S1_R2.reads.to.remove.txt
#
#   This wrapper (and the tool it wraps) required that both reads.to.remove.txt files be
#       applied to each fastq.  For example, to run remove.fastqs(), three inputs were required:
#           - fastq (e.g. S1_R1.fastq)
#           - forward reads.to.remove (e.g. S1_R1.reads.to.remove.txt)  
#           - reverse reads.to.remove (e.g. S1_R2.reads.to.remove.txt)  
#
#   This wrapper would then call (per the example above): 
#
#       remove.fastqs(S1_R1.fastq, S1_R1.reads.to.remove.txt, S1_R2.reads.to.remove.txt);
#
#   For the case above, the reads were not merged, so merged.reads should be set to FALSE.
#
#   If PEAR is used upstream in the pipeline (one of the first steps), then the forward and
#       reverse reads are merged into a single fastq, which is then used for the remainder
#       of the pipeline.
#   remove.fastqs() expects two reads.to.remove.txt files.  To keep revisions simple we 
#       simply modify this wrapper script and feed remove.fastqs() the same reads.to.remove.txt
#       file as input.  Functionally this gives us the results we want, without having to 
#       make large revisions to the code

#	Load required libraries
library(stringr);

#   Set this to TRUE if using both PEAR (single fastq input) AND you're using only
#       ONE reads-to-remove file
merged.reads <- FALSE;

#	Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);
working.dir <- arguments[1];

#   for debugging purposes
#working.dir <- "/Users/leyshock/Desktop/TCRseq/results/november/10_november/remove_fastq_testing/";

#	Specify relevant details
#   usage on ExaCloud
source("/mnt/lustre1/CompBio/genomic_resources/tcrseq/tcr_sequencing_tools/process_spikes/remove.spikes.R");
#   usage on personal laptop
#source("/Users/leyshock/Desktop/TCRseq/tools/tcr_sequencing_tools/tcr_sequencing_tools/process_spikes/remove.spikes.R");
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


