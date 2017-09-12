#	RScript that calculates number of unique clonotypes, Shannon entropy, and clonality for multiple clonotype files
#
#   Input is a directory containing one or more files with MiXCR output format
#   
#   Load necessary libraries
library(entropy);
library(data.table);
library(tcR)

#	Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);
#   Directory should contain all and only .txt files from exportClones, post-normalization with the following
#   format:
#     Clone count   Clone fraction    Clonal sequence(s)    AA. Seq. CDR3   Best V Hit    Best J Hit    V segments
#     J segments    Normalized clone count    Normalized clone fraction

clone.dir <- arguments[1];    # Typically .../dhaarini/DNAXXXXLC/normalization/normalized_clones/
count.dir <- arguments[2];    # Typically .../dhaarini/DNAXXXXLC/spike_counts/9bp/counts/
out.dir <- arguments[3];
batch_v <- arguments[4]


#	Examine the current directory for the files to process and sort them.
clone.files.in.dir <- list.files(clone.dir);
clone.files.in.dir <- clone.files.in.dir[order(as.numeric(gsub(".*_S|_.*$|\\..*$", '', clone.files.in.dir)))]
count.files.in.dir <- list.files(count.dir);
count.files.in.dir <- count.files.in.dir[order(as.numeric(gsub(".*_S|\\..*", '', count.files.in.dir)))]

if (length(clone.files.in.dir) != length(count.files.in.dir)) stop("Files do not match")


# Create empty arrays
calculated.entropies <- numeric(length(clone.files.in.dir));
unique.clones <- NULL
clonality <- NULL
max.clonal.freq <- numeric(length(count.files.in.dir));
norm.entropy <- numeric(length(clone.files.in.dir));
adaptive.clonality <- numeric(length(clone.files.in.dir));
adaptive.clonality <- NULL
max.clone.count <- numeric(length(clone.files.in.dir))
gini <- numeric(length(clone.files.in.dir))
true <- numeric(length(clone.files.in.dir))

for(i in 1:length(clone.files.in.dir))	{

    ###
    ### DATA
    ###
    #   get a clone file to process
    clone.curr.file <- clone.files.in.dir[i];
    clone.index <- strsplit(clone.curr.file, split = "_|\\.")[[1]][2]

    clone.curr.record <- fread(file.path(clone.dir, clone.curr.file),
                            check.names=FALSE,
                            stringsAsFactors=FALSE);
    #   get a count file to process
    count.curr.file <- count.files.in.dir[i];
    count.index <- strsplit(count.curr.file, split = "_|\\.")[[1]][2]

    count.curr.record <- fread(file.path(count.dir, count.curr.file))

    if (count.index != clone.index) stop("Mismatching files, check index: ", i)

    # Depending on if original or data_subset, we need a norm fraction column
    if ("New.norm.fraction" %in% colnames(clone.curr.record)){
      column <- "New.norm.fraction"
    } else {
      column <- "Normalized clone fraction"
    } # if

 #    column <- "Normalized.clone.fraction"

    ###
    ### CALCULATIONS
    ###

    # count number of lines in file, i.e. number of unique clonotypes
    unique.clones[i] <- clone.curr.record[,.N]

    #   calculate entropy
    calculated.entropies[i] <- entropy::entropy(clone.curr.record[[column]], method="ML", unit="log");


    #   calculate clonality
    clonality[i] <- 1 - (calculated.entropies[i] / log(unique.clones[i]))

    #   calculate clonality as the inverse of normalized shannon entropy
    	# Normalized shannon entropy
	norm.entropy[i] <- calculated.entropies[i] / log(unique.clones[i])
	# New clonality
    adaptive.clonality[i] <- 1 / norm.entropy[i]

    ## tcR Gini and True Diversity
    gini[i] <- gini(.data = clone.curr.record[[column]], .do.norm = F)
    true[i] <- diversity(.data = clone.curr.record[[column]], .do.norm = F)
	

    ## Change clone frequency column to a percentage
    clone.curr.record[[column]] <- clone.curr.record[[column]] * 100
    
    ##	Calculate Max. clonotype frequency
    max.clonal.freq[i] <- round(max(clone.curr.record[[column]]), digits = 4)

    #   Record maximum clone count
    max.clone.count[i] <- max(clone.curr.record$`Normalized clone count`)
#    max.clone.count[i] <- max(clone.curr.record$`Normalized.clone.count`)
    

    #   update progress
    if((i %%10) == 0)   {
        cat("Processing file ", i, " (out of ", length(count.files.in.dir), ")\n", sep="");
    }   #   fi

}	#	for i


#   create output data.frame
output.df <- data.frame(clone.files.in.dir, calculated.entropies, norm.entropy, unique.clones, clonality,
	     		adaptive.clonality, max.clonal.freq, max.clone.count, gini, true)

colnames(output.df) <- c("File", "Shannon Entropy", "Normalized Entropy", "Unique Clonotypes", "Clonality",
                         "Adaptive Clonality", "Max Clonal Freq", "Max Clone Count", "Gini_Index", "True_Diversity")

#   write output
#file.name <- "tcR.uniques.shannon.clonality.txt"
file.name <- paste0(batch_v, ".txt")
write.table(output.df, 
            file=paste(out.dir, file.name, sep = ''),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
