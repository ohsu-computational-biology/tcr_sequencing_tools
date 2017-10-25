### Analysis Script - Top Subset

### Calculate a variety of summary and diversity statistics for all of the normalized clonotype files
### Run the calculations on various subsets of the top clones
### Calculate Shannon Entropy, normalized Shannon Entropy, clonality, adaptive's clonality, unique clones, max clone count, gini index, true diversity

### Arguments
        ### clone.dir - directory containing normalized clone files
        ### count.dir - directory containing 9bp count files
        ### out.dir - directory to write analysis output
        ### old_v - TRUE - use old normalization method; FALSE - use new normalization method

### Load necessary libraries
library(entropy);
library(data.table);
library(tcR)

### Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);

clone.dir <- arguments[1];    # Typically .../dhaarini/DNAXXXXLC/normalization/normalized_clones/
count.dir <- arguments[2];    # Typically .../dhaarini/DNAXXXXLC/spike_counts/9bp/counts/
out.dir <- arguments[3];
old_v <- arguments[4]
divisions_v <- arguments[5] # How man divisions? comma-separated list of integers with no quotes or spaces (10,25,50,100)

### Examine the current directory for the files to process and sort them.
clone.files.in.dir <- list.files(clone.dir);
clone.files.in.dir <- clone.files.in.dir[order(as.numeric(gsub(".*_S|_.*$|\\..*$", '', clone.files.in.dir)))]

count.files.in.dir <- list.files(count.dir);
count.files.in.dir <- count.files.in.dir[order(as.numeric(gsub(".*_S|\\..*", '', count.files.in.dir)))]

if (length(clone.files.in.dir) != length(count.files.in.dir)) stop("Files do not match")

### Get divisions
divisions_v <- sapply(strsplit(divisions_v, split = ',')[[1]], function(x) as.numeric(x), USE.NAMES = F)

### Run analysis for each division
for (j in 1:length(divisions_v) {

    ## Get divisions
    currDiv_v <- divisions_v[j]

    ## Create arrays
    clonality <- numeric(length(clone.files.in.dir))
    calculated.entropies <- clonality
    unique.clones <- clonality
    max.clonal.freq <- clonality
    norm.entropy <- clonality
    adaptive.clonality <- clonality
    max.clone.count <- clonality
    gini <- clonality
    trueD <- clonality

    ## Update
    print(c("Working on division: ", currDiv_v))

    for(i in 1:length(clone.files.in.dir))	{

        ###
        ### DATA
        ###
    
        ## Get a clone file to process
        clone.curr.file <- clone.files.in.dir[i];
        clone.index <- strsplit(clone.curr.file, split = "_|\\.")[[1]][2]

        clone.curr.record <- fread(file.path(clone.dir, clone.curr.file))

        ## Get a count file to process
        count.curr.file <- count.files.in.dir[i];
        count.index <- strsplit(count.curr.file, split = "_|\\.")[[1]][2]

        count.curr.record <- fread(file.path(count.dir, count.curr.file))

        if (count.index != clone.index) stop("Mismatching files, check index: ", i)

        ## Get fraction column
        if (old_v) {
            column_v <- "Normalized clone fraction"
            count_v <- "Normalized clone count"
        } else {
            column_v <- "nb.clone.fraction"
            count_v <- "nb.clone.count"
        } # fi

        ## Sort and subset
	clone.curr.record <- clone.cur.record[order(clone.curr.record[[column_v]], decreasing = T),]
        clone.curr.record <- clone.curr.record[1:currDiv_v,]

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
        print(paste0(clone.index, " Gini and true"))
        gini[i] <- gini(.data = clone.curr.record[[column_v]], .do.norm = T)
        trueD[i] <- diversity(.data = clone.curr.record[[column_v]], .do.norm = T)

        ## Change clone frequency column to a percentage
        clone.curr.record[[column_v]] <- clone.curr.record[[column_v]] * 100
    
        ## Calculate Max. clonotype frequency
        max.clonal.freq[i] <- round(max(clone.curr.record[[column_v]]), digits = 4)

        ## Record maximum clone count
        max.clone.count[i] <- max(clone.curr.record[[count_v]])
    
        ## Update progress
        if((i %%10) == 0)   {
            cat("Processing file ", i, " (out of ", length(count.files.in.dir), ")\n", sep="");
        } # fi

    } # for i

    ## Create output data.frame
    output.df <- data.frame(clone.files.in.dir, calculated.entropies, norm.entropy, unique.clones, clonality,
	     		    adaptive.clonality, max.clonal.freq, max.clone.count, gini, trueD)

    colnames(output.df) <- c("File", "Shannon Entropy", "Normalized Entropy", "Unique Clonotypes", "Clonality",
                             "Adaptive Clonality", "Max Clonal Freq", "Max Clone Count", "Gini_Index", "True_Diversity")

    ## Write output
    file.name <- paste0(batch_v, "_", currDiv_v, ".txt")
    write.table(output.df, 
                file=paste(out.dir, file.name, sep = ''),
                quote=FALSE,
                sep="\t",
                row.names=FALSE)
