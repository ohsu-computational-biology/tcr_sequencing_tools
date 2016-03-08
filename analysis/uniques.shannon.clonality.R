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
count.dir <- arguments[2];    # Typically .../dhaarini/DNAXXXXLC/spike_counts/9bp/counts/


#	Examine the current directory for the files to process and sort them.
clone.files.in.dir <- list.files(clone.dir);
clone.files.in.dir <- clone.files.in.dir[order(as.numeric(gsub(".*_S|_alignment_.*", '', clone.files.in.dir)))]
count.files.in.dir <- list.files(count.dir);
count.files.in.dir <- count.files.in.dir[order(as.numeric(gsub(".*_S|//..*", '', count.files.in.dir)))]


# Create empty arrays
calculated.entropies <- numeric(length(clone.files.in.dir));
unique.clones <- NULL
clonality <- NULL
spike.counts <- numeric(length(count.files.in.dir));
spike.percent <- numeric(length(count.files.in.dir));
max.clonal.freq <- numeric(length(count.files.in.dir));

for(i in 1:length(clone.files.in.dir))	{
    #   get a clone file to process
    clone.curr.file <- clone.files.in.dir[i];

    clone.curr.record <- read.delim(file.path(clone.dir, clone.curr.file),
                            check.names=FALSE,
                            stringsAsFactors=FALSE);
   #   get a count file to process
   count.curr.file <- count.files.in.dir[i];

   count.curr.record <- read.delim(file.path(count.dir, count.curr.file), sep = ',',
   		     		check.names=FALSE, stringsAsFactors=FALSE);

    # count number of lines in file, i.e. number of unique clonotypes
    curr.clone.system.call <- paste("wc -l ", 
                                  clone.dir, clone.files.in.dir[i], 
                                  " | awk '{print $1}'",
                                  sep="");
    unique.clones[i] <- as.numeric(system(curr.clone.system.call, intern=TRUE)) - 1;

    #   calculate entropy
    calculated.entropies[i] <- entropy(clone.curr.record$"Normalized clone fraction", method="ML", unit="log");

    #   calculate clonality
    clonality[i] <- 1 - (calculated.entropies[i] / log(unique.clones[i]))

    #   Record Spike Counts
    spike.counts[i] <- count.curr.record[1,5];

    #   Record % spiked reads
    spike.percent[i] <- spike.counts[i] / sum(clone.curr.record$"Normalized clone count") * 100

    #	Calculate Max. clonotype frequency
    max.clonal.freq[i] <- max(clone.curr.record$`Normalized clone fraction`) * 100

    #   update progress
    if((i %%10) == 0)   {
        cat("Processing file ", i, " (out of ", length(files.in.dir), ")\n", sep="");
    }   #   fi

}	#	for i

#   create output data.frame
output.df <- data.frame(clone.files.in.dir, calculated.entropies, unique.clones, clonality,
	     		spike.counts, spike.percent, max.clonal.freq);

colnames(output.df) <- c("File", "Shannon Entropy", "Unique Clonotypes", "Clonality", "Spiked Reads",
		       	"Percent Spiked Reads", "Max Clonal Freq");

#   write output
write.table(output.df, 
            file="uniques.shannon.clonality.txt",
            quote=FALSE,
            sep=",",
            row.names=FALSE)
