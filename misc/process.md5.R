#   This script:
#		1.  calculates the MD5 sums for all files in a directory
#		2.  sorts the md5sum file provided by the Core
#   
#   Input:  directory containing one or more .fastq.gz files
#   Ouputs:  tab-separated files with the following format:
#
#           md5sum  file.name.fastq.gz
#
#   Run this script from the command line as following:
#
#           ~% Rscript compute.md5.Rscript "input_directory/"
#
#	If the files have been transferred without error, diffing the files as
#		specified in the program output should return an empty result

#   Load appropriate libraries
.libPaths("/home/exacloud/lustre1/CompBio/lib/R/library")
library(digest);

#	get appropriate command line directories
arguments <- commandArgs(trailingOnly=TRUE);
input.directory <- arguments[1];

#   get the files to be processed, and specify the output paths
fastq.files.to.process <- list.files(input.directory, pattern="*.fastq.gz");
pathed.fastq.files.to.process <- paste(input.directory, fastq.files.to.process, sep="");
pathed.digest.file.to.process <- paste(input.directory, "md5sum.txt", sep="");
output.calculated.digest.file.location <- paste(input.directory, "calculated.md5.sums.txt", sep=""); 
output.sorted.digest.file.location <- paste(input.directory, "md5sum.sorted.txt", sep=""); 

##################################################################################
#	Calculate md5 sums for the provided fastq.gz files
##################################################################################

#   create object to store results 
output.digests <- character(length(fastq.files.to.process));

cat("Computing MD5 sums for ", length(fastq.files.to.process), " fastq.gz files\n", sep="");

for(i in 1:length(fastq.files.to.process))    {
    output.digests[i] <- digest(algo="md5",
                                serialize=FALSE,
                                file=pathed.fastq.files.to.process[i]);
	if((i %% 10) == 0)	{
		cat("Digesting file number ", i, "\n", sep="");
	}	#	fi
}   #   for i


#   create output data.frame
output.df <- data.frame(output.digests, fastq.files.to.process);

cat("Writing calculated digest file to: ", output.calculated.digest.file.location, "\n", sep="");

#   output results
write.table(output.df, 
            file=output.calculated.digest.file.location,
            quote=FALSE,
            sep=" ",
            row.names=FALSE,
            col.names=FALSE);


##################################################################################
#	Sort the md5 sums provided by the core
##################################################################################


provided.digests <- read.csv(pathed.digest.file.to.process,
							header=FALSE,
							sep=" ",
							stringsAsFactors=FALSE);

provided.digests$V2 <- NULL;

provided.digests <- provided.digests[order(provided.digests$V3),];

sorted.digest.location <- paste(input.directory, "md5sum.sorted.txt", sep=""); 
cat("Writing sorted digest file to: ", output.sorted.digest.file.location, "\n", sep="");

#   output results
write.table(provided.digests, 
            file=output.sorted.digest.file.location,
            quote=FALSE,
            sep=" ",
            row.names=FALSE,
            col.names=FALSE);

cat("Diff ", output.calculated.digest.file.location, " and ", output.sorted.digest.file.location, "\n", sep=""); 
