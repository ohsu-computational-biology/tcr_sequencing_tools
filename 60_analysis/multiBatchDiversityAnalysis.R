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

### Get all of the files
clone.files.in.dir <- list.files(clone.dir);

### Get batches
batches_v <- unique(sapply(clone.files.in.dir, function(x) strsplit(x, split = "_")[[1]][1]))

### Sort files by batch and then by sample within batch
clone.files.in.dir <- unname(unlist(sapply(batches_v, function(x) {
	y <- clone.files.in.dir[grep(x, clone.files.in.dir)]
	z <- y[order(as.numeric(gsub(".*_S|_.*$|\\..*$}", '', y)))]
	return(z)})))


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
batchCol_v <- character(length(clone.files.in.dir))

for(i in 1:length(clone.files.in.dir))	{

    ###
    ### DATA
    ###

    ## Get a clone file to process
    clone.curr.file <- clone.files.in.dir[i];
    clone.index <- grep("S[0-9]+", strsplit(clone.curr.file, split = "_|\\.")[[1]], value = T)
    clone.curr.record <- fread(file.path(clone.dir, clone.curr.file))

    ## Get batch
    currBatch_v <- strsplit(clone.curr.file, split = "_")[[1]][1]
    batchCol_v[i] <- currBatch_v

    ## Get columns
    column_v <- "normFreq"
    count_v <- "normC"

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
file.name <- paste0("multiBatch_", paste(batches_v, collapse = "_"), "_divAnalysis.txt")
write.table(output.df, 
            file=file.path(out.dir, file.name),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
