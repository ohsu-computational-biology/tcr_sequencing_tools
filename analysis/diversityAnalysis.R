### Analysis Script

### Calculate a variety of summary and diversity statistics for all of the normalized clonotype files

### Calculate Shannon Entropy, normalized Shannon Entropy, clonality, adaptive's clonality, unique clones, max clone count, gini index, true diversity

### Arguments
	### clone.dir - directory containing normalized clone files
	### out.dir - directory to write analysis output
        ### old_v - TRUE - use old normalization method; FALSE - use new normalization method
###   Load necessary libraries
suppressMessages(library(entropy))
suppressMessages(library(data.table))
suppressMessages(library(tcR))

### Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);

clone.dir <- arguments[1];    # Typically .../dhaarini/DNAXXXXLC/normalization/normalized_clones/
out.dir <- arguments[2];
old_v <- arguments[3]

### Examine the current directory for the files to process and sort them.
clone.files.in.dir <- list.files(clone.dir);
clone.files.in.dir <- clone.files.in.dir[order(as.numeric(gsub(".*_S|_.*$|\\..*$", '', clone.files.in.dir)))]

### Get batch
batch_v <- strsplit(clone.files.in.dir[1], split = "_")[[1]][1]

### Create empty arrays
calculated.entropies <- numeric(length(clone.files.in.dir));
unique.clones <- numeric(length(clone.files.in.dir))
clonality <- numeric(length(clone.files.in.dir))
max.clonal.freq <- numeric(length(clone.files.in.dir));
norm.entropy <- numeric(length(clone.files.in.dir));
adaptive.clonality <- numeric(length(clone.files.in.dir));
adaptive.clonality <- numeric(length(clone.files.in.dir))
max.clone.count <- numeric(length(clone.files.in.dir))
gini <- numeric(length(clone.files.in.dir))
true <- numeric(length(clone.files.in.dir))

for(i in 1:length(clone.files.in.dir))	{

    ###
    ### DATA
    ###

    ## Get a clone file to process
    clone.curr.file <- clone.files.in.dir[i];
    clone.index <- grep("S[0-9]+", strsplit(clone.curr.file, split = "_|\\.")[[1]], value = T)
    clone.curr.record <- fread(file.path(clone.dir, clone.curr.file))

    ## Get fraction column
    if (old_v) {
        column_v <- "Normalized clone fraction"
        count_v <- "Normalized clone count"
    } else {
	column_v <- "nb.clone.fraction"
        count_v <- "nb.clone.count"
    } # fi

    ###
    ### CALCULATIONS
    ###

    ## Count number of lines in file, i.e. number of unique clonotypes
    unique.clones[i] <- clone.curr.record[,.N]

    ## Calculate entropy
    calculated.entropies[i] <- entropy::entropy(clone.curr.record[[column_v]], method="ML", unit="log");

    ## Calculate clonality
    clonality[i] <- 1 - (calculated.entropies[i] / log(unique.clones[i]))

    ## Calculate clonality as the inverse of normalized shannon entropy
    ## Normalized shannon entropy
    norm.entropy[i] <- calculated.entropies[i] / log(unique.clones[i])
    ## New clonality
    adaptive.clonality[i] <- 1 / norm.entropy[i]

    ## tcR Gini and True Diversity
    gini[i] <- gini(.data = clone.curr.record[[column_v]], .do.norm = F)
    true[i] <- diversity(.data = clone.curr.record[[column_v]], .do.norm = F)
	
    ## Change clone frequency column to a percentage
    clone.curr.record[[column_v]] <- clone.curr.record[[column_v]] * 100
    
    ## Calculate Max. clonotype frequency
    max.clonal.freq[i] <- round(max(clone.curr.record[[column_v]]), digits = 4)

    ## Record maximum clone count
    max.clone.count[i] <- max(clone.curr.record[[count_v]])

    ##   update progress
    if((i %%10) == 0)   {
        cat("Processing file ", i, " (out of ", length(clone.files.in.dir), ")\n", sep="");
    } # fi

} # for i

### Create output data.frame
output.df <- data.frame(clone.files.in.dir, calculated.entropies, norm.entropy, unique.clones, clonality,
	     		adaptive.clonality, max.clonal.freq, max.clone.count, gini, true)

colnames(output.df) <- c("File", "Shannon Entropy", "Normalized Entropy", "Unique Clonotypes", "Clonality",
                         "Adaptive Clonality", "Max Clonal Freq", "Max Clone Count", "Gini_Index", "True_Diversity")

### Write output
file.name <- paste0(batch_v, ".txt")
write.table(output.df, 
            file=paste(out.dir, file.name, sep = ''),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
