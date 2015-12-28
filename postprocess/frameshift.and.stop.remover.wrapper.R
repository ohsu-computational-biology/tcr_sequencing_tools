
#
#	Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);
working.dir <- arguments[1];
remove.records <- arguments[2];

#	Specify path to script
#	For use on ExaCloud
source("/mnt/lustre1/CompBio/genomic_resources/tcrseq/tcr_sequencing_tools/postprocess/frameshift.and.stop.remover.R");
#	For use on laptop
#source("/Users/leyshock/Desktop/TCRseq/tools/tcr_sequencing_tools/postprocess/frameshift.and.stop.remover.R");

#	Examine the current directory for the files to process
files.in.dir <- list.files(working.dir);

for(i in 1:length(files.in.dir))  {
	#	call function
	postprocess.clone.table(file.path(working.dir, files.in.dir[i]), remove.records);
    cat("Processed: ", i, " (of ", length(files.in.dir), ")\n",sep="");

}	#	for i

