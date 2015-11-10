
#
#	Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);
working.dir <- arguments[1];

#	Specify relevant details
##source("/mnt/lustre1/CompBio/genomic_resources/tcrseq/optimized.remove.spikes.R");
source("/Users/leyshock/Desktop/TCRseq/tools/tcr_sequencing_tools/frameshift.and.stop.remover.R");

#	Examine the current directory for the files to process
files.in.dir <- list.files(working.dir);

for(i in 1:length(files.in.dir))  {
	#	call function
	remove.records(files.in.dir[i]);
    cat("Processed: ", i, " (of ", length(files.in.dir), ")\n",sep="");

}	#	for i

