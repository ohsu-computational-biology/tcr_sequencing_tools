#	RScript that calculates number of unique clonotypes, Shannon entropy, and clonality for multiple clonotype files
#
#   Input is a directory containing one or more files with MiXCR output format
#   
#   Load necessary libraries
library(entropy);

#	Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);
#   Directory should contain all and only .txt files from exportClones, post-normalization with the following
#   format:
#     Clone count   Clone fraction    Clonal sequence(s)    AA. Seq. CDR3   Best V Hit    Best J Hit    V segments
#     J segments    Normalized clone count    Normalized clone fraction

clone.dir <- arguments[1];    # Typically .../dhaarini/DNAXXXXLC/normalization/normalized_clones/

#	Examine the current directory for the files to process
clone.files.in.dir <- list.files(clone.dir);

calculated.entropies <- numeric(length(clone.files.in.dir));
unique.clones <- NULL
clonality <- NULL

for(i in 1:length(clone.files.in.dir))	{
    #   get a clone file to process
    clone.curr.file <- clone.files.in.dir[i];

    clone.curr.record <- read.delim(file.path(clone.dir, clone.curr.file),
                            check.names=FALSE,
                            stringsAsFactors=FALSE);
    # count number of lines in file, i.e. number of unique clonotypes
    curr.clone.system.call <- paste("wc -l ", 
                                  clone.dir, clone.files.in.dir[i], 
                                  " | awk '{print $1}'",
                                  sep="");
    unique.clones[i] <- as.numeric(system(curr.clone.system.call, intern=TRUE)) - 1;

    #   calculate entropy
    calculated.entropies[i] <- entropy(clone.curr.record$"Clone fraction", method="ML", unit="log");

    #   calculate clonality
    clonality[i] <- 1 - calculated.entropies[i] / log(unique.clones[i])

    #   update progress
    if((i %%10) == 0)   {
        cat("Processing file ", i, " (out of ", length(clone.files.in.dir), ")\n", sep="");
    }   #   fi

}	#	for i

#   create output data.frame
output.df <- data.frame(clone.files.in.dir, calculated.entropies, unique.clones, clonality);

#   write output
write.table(output.df, 
            file="uniques.shannon.clonality.txt",
            quote=FALSE,
            sep=",",
            row.names=FALSE)
