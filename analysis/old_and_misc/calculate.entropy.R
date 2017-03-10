#	RScript that calculates Shannon entropy for multiple clonotype files
#
#   Input is a directory containing one or more files with MiXCR output format
#   
#   Load necessary libraries
.libPaths("/mnt/lustre1/CompBio/lib/R/Library")
library(entropy); # requires package "entropy" as well

#	Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);
#   Directory should contain all and only QC files for a given "batch"
working.dir <- arguments[1];

#	Examine the current directory for the files to process
files.in.dir <- list.files(working.dir);

calculated.entropies <- numeric(length(files.in.dir));

for(i in 1:length(files.in.dir))	{
    #   get a QC file to process
    curr.file <- files.in.dir[i];

    curr.record <- read.delim(file.path(working.dir, curr.file),
                            check.names=FALSE,
                            stringsAsFactors=FALSE);

    #   calculate entropy
    #   Note:  The library's code throw's warnings about the input not summing to 1, even
    #       though we're using the .do.norm parameter in the function call
    calculated.entropies[i] <- entropy(curr.record$"Normalized clone fraction", method="ML", unit="log");

    #   update progress
    if((i %%10) == 0)   {
        cat("Processing file ", i, " (out of ", length(files.in.dir), ")\n", sep="");
    }   #   fi

}	#	for i

#   create output data.frame
output.df <- data.frame(files.in.dir, calculated.entropies);

write.table(output.df, 
            file="calculated.entropies.txt",
            quote=FALSE,
            sep=",",
            row.names=FALSE)


