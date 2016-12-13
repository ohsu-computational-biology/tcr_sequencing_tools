#	RScript that calculates number of unique clonotypes, Shannon entropy, and clonality for multiple clonotype files
#
#   Input is a directory containing one or more files with MiXCR output format
#   
#   Load necessary libraries
.libPaths("/home/exacloud/lustre1/CompBio/lib/R/library")
library(entropy);

#	Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);
#   Directory should contain all and only .txt files from exportClones, post-normalization with the following
#   format:
#     Clone count   Clone fraction    Clonal sequence(s)    AA. Seq. CDR3   Best V Hit    Best J Hit    V segments
#     J segments    Normalized clone count    Normalized clone fraction

clone.dir <- arguments[1];    # Typically .../dhaarini/DNAXXXXLC/normalization/normalized_clones/
#count.dir <- arguments[2];    # Typically .../dhaarini/DNAXXXXLC/spike_counts/9bp/counts/
batch <- arguments[2]
out.dir <- arguments[3];


#	Examine the current directory for the files to process and sort them.
clone.files.in.dir <- list.files(clone.dir);
if (batch == "RNA140128LC") {
   clone.files.in.dir <- clone.files.in.dir[order(as.numeric(gsub(".*LC_|_exported.*", '', clone.files.in.dir)))]
   } else {
   clone.files.in.dir <- clone.files.in.dir[order(as.numeric(gsub(".*LC_[A-Z]|_exported.*", '', clone.files.in.dir)))]
} # fi
#clone.files.in.dir <- clone.files.in.dir[order(as.numeric(gsub(".*_S|_exported.*", '', clone.files.in.dir)))]
#count.files.in.dir <- list.files(count.dir);
#count.files.in.dir <- count.files.in.dir[order(as.numeric(gsub(".*_S|\\..*", '', count.files.in.dir)))]


# Create empty arrays
calculated.entropies <- numeric(length(clone.files.in.dir));
unique.clones <- NULL
clonality <- NULL
#spike.counts <- numeric(length(count.files.in.dir));
#spike.percent <- numeric(length(count.files.in.dir));
max.clonal.freq <- numeric(length(clone.files.in.dir));
norm.entropy <- numeric(length(clone.files.in.dir));
adaptive.clonality <- numeric(length(clone.files.in.dir));
adaptive.clonality <- NULL
max.clone.count <- numeric(length(clone.files.in.dir))
top.10 <- data.frame(matrix(nrow = 10, ncol = length(clone.files.in.dir)))
top.25 <- data.frame(matrix(nrow = 25, ncol = length(clone.files.in.dir)))
top.50 <- data.frame(matrix(nrow = 50, ncol = length(clone.files.in.dir)))

for(i in 1:length(clone.files.in.dir))	{

    ###
    ### DATA
    ###
    #   get a clone file to process
    clone.curr.file <- clone.files.in.dir[i];

    # Check size
    clone.curr.size <- file.size(paste(clone.dir, clone.curr.file, sep = ''))

    if (clone.curr.size == 0){next}

    clone.curr.record <- read.delim(file.path(clone.dir, clone.curr.file),
                            check.names=FALSE,
                            stringsAsFactors=FALSE);
    #   get a count file to process
#    count.curr.file <- count.files.in.dir[i];

#    count.curr.record <- read.delim(file.path(count.dir, count.curr.file), sep = ',',
#   		     		check.names=FALSE, stringsAsFactors=FALSE);

    # Depending on if original or data_subset, we need a norm fraction column
#    if ("New.norm.fraction" %in% colnames(clone.curr.record)){
#      column <- "New.norm.fraction"
#    } else {
#      column <- "Normalized clone fraction"
#    } # if

 #    column <- "Normalized.clone.fraction"
      column <- "Clone fraction"

    ###
    ### CALCULATIONS
    ###

    # count number of lines in file, i.e. number of unique clonotypes
    unique.clones[i] <- length(clone.curr.record[,1])

    #   calculate entropy
    calculated.entropies[i] <- entropy(clone.curr.record[[column]], method="ML", unit="log");


    #   calculate clonality
    clonality[i] <- 1 - (calculated.entropies[i] / log(unique.clones[i]))

    #   calculate clonality as the inverse of normalized shannon entropy
    	# Normalized shannon entropy
	norm.entropy[i] <- calculated.entropies[i] / log(unique.clones[i])
	# New clonality
	adaptive.clonality[i] <- 1 / norm.entropy[i]
	

    #   Record Spike Counts
#    spike.counts[i] <- count.curr.record[1,5];

    #   Record % spiked reads
#    spike.percent[i] <- spike.counts[i] / sum(clone.curr.record$"Clone count") * 100
#    spike.percent[i] <- spike.counts[i] / sum(clone.curr.record$"Normalized clone count") * 100

    #	Calculate Max. clonotype frequency
    max.clonal.freq[i] <- max(clone.curr.record[[column]]) * 100

    #   Record maximum clone count
#    max.clone.count[i] <- max(clone.curr.record$`Normalized clone count`)
    max.clone.count[i] <- max(clone.curr.record$`Clone count`)

   #  Record frequencies for top 10 and top 25 clones
    clone.curr.record <- clone.curr.record[order(clone.curr.record[[column]], decreasing = T),]
    top.10[,i] <- clone.curr.record[[column]][1:10]
    top.25[,i] <- clone.curr.record[[column]][1:25]
    top.50[,i] <- clone.curr.record[[column]][1:50]

    #   update progress
    if((i %%10) == 0)   {
        cat("Processing file ", i, " (out of ", length(clone.files.in.dir), ")\n", sep="");
    }   #   fi

}	#	for i

# Summarize top10 and top25 data
top.10.summary <- t(apply(top.10, 2, function(x) c(mean(x), median(x), sum(x))))
top.25.summary <- t(apply(top.25, 2, function(x) c(mean(x), median(x), sum(x))))
top.50.summary <- t(apply(top.50, 2, function(x) c(mean(x), median(x), sum(x))))

#   create output data.frame
output.df <- data.frame(clone.files.in.dir, calculated.entropies, norm.entropy, unique.clones, clonality,
	     		adaptive.clonality, max.clonal.freq, max.clone.count,
			top.10.summary, top.25.summary, top.50.summary);

colnames(output.df) <- c("File", "Shannon Entropy", "Normalized Entropy", "Unique Clonotypes", "Clonality",
		        "Adaptive Clonality", "Max Clonal Freq",
			"Max Clone Count", "top10.mean", "top10.median", "top10.sum", "top25.mean", "top25.median",
			"top25.sum", "top50.mean", "top50.median", "top50.sum");

#   write output
file.name <- "uniques.shannon.clonality.txt"
write.table(output.df, 
            file=paste(out.dir, file.name, sep = ''),
            quote=FALSE,
            sep=",",
            row.names=FALSE)
