library(stringr);

tester <- function(working.dir) {
#
#	Get command-line arguments
#arguments <- commandArgs(trailingOnly=TRUE);
#working.dir <- arguments[1];

#	Specify relevant details
##source("/mnt/lustre1/CompBio/genomic_resources/tcrseq/optimized.remove.spikes.R");
source("/Users/leyshock/Desktop/TCRseq/tools/tcr_sequencing_tools/normalize.clonotype.counts.R");
counts.suffix <- "_counts.txt";
clones.suffix <- "_exported_clones.txt";

#	Examine the current directory for the files to process
files.in.dir <- list.files(working.dir);

#	Prepare a list of files to process
clone.files <- str_detect(files.in.dir, clones.suffix);
clone.files <- files.in.dir[clone.files];
clone.files <- paste(working.dir, clone.files, sep="");
count.files <- str_detect(files.in.dir, counts.suffix);
count.files <- files.in.dir[count.files];
count.files <- paste(working.dir, count.files, sep="");
file.roots.1 <- unique(clone.files);
file.roots.2 <- unique(count.files);

#	error-check
if(length(file.roots.1) != length(file.roots.2))
	stop("Mismatch between number of count.files and number of clone.files in directory")
for(i in 1:length(count.files))	{

	#	call function
	normalize.clonotype.counts(clone.files[i], count.files[i], i);

	if((i %% 10) == 0)	{
		cat("Processed: ", i, " (of ", length(count.files), ")\n",sep="");
	}	#	fi

}	#	for i

}
