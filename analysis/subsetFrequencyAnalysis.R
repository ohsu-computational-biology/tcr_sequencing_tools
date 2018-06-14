### Analysis Script - Top Subset

### Calculate a variety of summary and diversity statistics for all of the normalized clonotype files
### Run the calculations on various subsets of the top clones
### NOTE - this is not 'frequency groups' in terms of Small, Medium, Large, Hyper, but just different subsets of the top clones
### Calculate Shannon Entropy, normalized Shannon Entropy, clonality, adaptive's clonality, unique clones, max clone count, gini index, true diversity

### Arguments
        ### clone.dir - directory containing normalized clone files
        ### out.dir - directory to write analysis output
        ### old_v - TRUE - use old normalization method; FALSE - use new normalization method

### Load necessary libraries
suppressMessages(library(entropy))
suppressMessages(library(data.table))

### Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);

clone.dir <- arguments[1];    # Typically .../dhaarini/DNAXXXXLC/normalization/normalized_clones/
out.dir <- arguments[2];
old_v <- arguments[3]
divisions_v <- arguments[4] # How man divisions? comma-separated list of integers with no quotes or spaces (10,25,50,100)

### Examine the current directory for the files to process and sort them.
clone.files.in.dir <- list.files(clone.dir);
clone.files.in.dir <- clone.files.in.dir[order(as.numeric(gsub(".*_S|_.*$|\\..*$", '', clone.files.in.dir)))]

### Get batch
batch_v <- strsplit(clone.files.in.dir[1], split = "_")[[1]][1]

### Get divisions
divisions_v <- sapply(strsplit(divisions_v, split = ',')[[1]], function(x) as.numeric(x), USE.NAMES = F)

### Empty output table
totalFreq_dt <- NULL

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

        ## Sort
	clone.curr.record <- clone.curr.record[order(clone.curr.record[[column_v]], decreasing = T),]

	## Change clone frequency into a percentage
	clone.curr.record[[column_v]] <- clone.curr.record[[column_v]] * 100

	## Emtpy variables
	totalOut_v <- NULL
	totalNames_v <- NULL

	## For each division, calculate the mean and sum of the clonal frequencies
	for (j in 1:length(divisions_v)){
	    ## Get division
	    currDiv_v <- divisions_v[j]

	    ## Subset and also remove NAs if the sample has fewere clones than the division
	    currData_dt <- clone.curr.record[1:currDiv_v,][!is.na(`V segments`)]
	    numClones_v <- nrow(currData_dt)

	    ## Warn about too small
	    if (numClones_v < currDiv_v){
		print(sprintf("Sample %s only has %d clones, using that instead of top %d clones", clone.index, numClones_v, currDiv_v))
	    } # fi

	    ## Get mean and sum of frequencies
	    currMeanFreq_v <- mean(currData_dt[[column_v]])
	    currSumFreq_v <- sum(currData_dt[[column_v]])

	    ## Combine together
            currNames_v <- c(paste0(currDiv_v, "_clonesUsed"), paste0(currDiv_v, "_meanFreq"), paste0(currDiv_v, "_sumFreq"))
	    currOut_v <- c(nrow(currData_dt), currMeanFreq_v, currSumFreq_v)

	    ## Add to total output
	    totalOut_v <- c(totalOut_v, currOut_v)
	    totalNames_v <- c(totalNames_v, currNames_v)
	} # for j

	## Combine with total
	totalFreq_dt <- rbind(totalFreq_dt, totalOut_v)

        ## Update progress
        if((i %%10) == 0)   {
            cat("Processing file ", i, " (out of ", length(clone.files.in.dir), ")\n", sep="");
        } # fi

    } # for i

## Create final output
output.df <- cbind(clone.files.in.dir, totalFreq_dt)
colnames(output.df) <- c("File", totalNames_v)

## Write output
file.name <- paste0(batch_v, "_topCloneFreqs.txt")
write.table(output.df, 
            file=paste(out.dir, file.name, sep = ''),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
