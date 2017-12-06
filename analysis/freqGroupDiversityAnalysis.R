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
suppressMessages(library(entropy));
suppressMessages(library(data.table));
suppressMessages(library(tcR))

### Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);

clone.file <- arguments[1];    # Typically .../dhaarini/DNAXXXXLC/freqGroups/BATCH_full_clones.txt/
out.dir <- arguments[2];
old_v <- arguments[3]
divisions_v <- arguments[4] # Which frequency groups to divide (comma-sep, no spaces)

### Read data
cloneData_dt <- fread(clone.file)

### Get divisions
### If using a specific division file, the divisions_v argument will be NA, and we'll want to grab the division from the file name
### If using the "full" file, one or many divisions may be specified, spearated by commas
if (is.na(divisions_v)) {
    divisions_v <- strsplit(clone.file, split = "_")[[1]][2]
} else {
    divisions_v <- unlist(strsplit(divisions_v, split = ','))
} # fi

### Get samples
sampCol_v <- grep("ample", colnames(cloneData_dt), value = T)
samples_v <- unique(cloneData_dt[[sampCol_v]])

### Get batch
batch_v <- strsplit(basename(clone.file), split = "_")[[1]][1]

### Run analysis for each division
for (j in 1:length(divisions_v)) {

    ## Get divisions
    currDiv_v <- divisions_v[j]

    ## Subset Data
    currDivData_dt <- cloneData_dt[Div == currDiv_v,]
    
    ## Update
    print(c("Working on division: ", currDiv_v))
    print(c("Number of clones in this division: ", nrow(currDivData_dt)))

    ## Update samples
    print(c("Original samples: ", samples_v))
    subSamples_v <- unique(currDivData_dt[[sampCol_v]])
    print(c("Samples in this division: ", subSamples_v))

    ## Create arrays
    clonality <- numeric(length(subSamples_v))
    calculated.entropies <- clonality
    unique.clones <- clonality
    max.clonal.freq <- clonality
    norm.entropy <- clonality
    adaptive.clonality <- clonality
    max.clone.count <- clonality
    gini <- clonality
    trueD <- clonality

    ## Get fraction column
    if (old_v) {
        column_v <- "Normalized clone fraction"
        count_v <- "Normalized clone count"
    } else {
        column_v <- "nb.clone.fraction"
        count_v <- "nb.clone.count"
    } # fi

    for(i in 1:length(subSamples_v)) {

        ## Get current sample
        currSample_v <- subSamples_v[i]
        currData_dt <- currDivData_dt[get(sampCol_v) == currSample_v,]

        ## Count number of lines in file, i.e. number of unique clonotypes
        unique.clones[i] <- currData_dt[,.N]

        ## Calculate entropy
        calculated.entropies[i] <- entropy::entropy(currData_dt[[column_v]], method="ML", unit="log");

        ## Calculate clonality
        clonality[i] <- 1 - (calculated.entropies[i] / log(unique.clones[i]))

        ## Calculate clonality as the inverse of normalized shannon entropy
    	## Normalized shannon entropy
	norm.entropy[i] <- calculated.entropies[i] / log(unique.clones[i])
	## New clonality
        adaptive.clonality[i] <- 1 / norm.entropy[i]

        ## tcR Gini and True Diversity
        gini[i] <- gini(.data = currData_dt[[column_v]], .do.norm = T)
        trueD[i] <- diversity(.data = currData_dt[[column_v]], .do.norm = T)

        ## Change clone frequency column to a percentage
        currData_dt[[column_v]] <- currData_dt[[column_v]] * 100
    
        ## Calculate Max. clonotype frequency
        max.clonal.freq[i] <- round(max(currData_dt[[column_v]]), digits = 4)

        ## Record maximum clone count
        max.clone.count[i] <- max(currData_dt[[count_v]])
    
        ## Update progress
        if((i %%10) == 0)   {
            cat("Processing file ", i, " (out of ", length(subSamples_v), ")\n", sep="");
        } # fi

    } # for i

    ## Create output data.frame
    output.df <- data.frame(subSamples_v, calculated.entropies, norm.entropy, unique.clones, clonality,
	     		    adaptive.clonality, max.clonal.freq, max.clone.count, gini, trueD)

    colnames(output.df) <- c("Sample", "Shannon Entropy", "Normalized Entropy", "Unique Clonotypes", "Clonality",
                             "Adaptive Clonality", "Max Clonal Freq", "Max Clone Count", "Gini_Index", "True_Diversity")

    ## Write output
    file.name <- paste0(batch_v, "_", currDiv_v, ".txt")
    write.table(output.df, 
                file=paste(out.dir, file.name, sep = ''),
                quote=FALSE,
                sep="\t",
                row.names=FALSE)
} # for j
