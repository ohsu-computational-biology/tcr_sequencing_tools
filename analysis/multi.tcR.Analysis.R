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
divisions_v <- arguments[4];  # True means to get sum/mean/median of top clones


#	Examine the current directory for the files to process and sort them.
clone.files.in.dir <- list.files(clone.dir);
clone.files.in.dir <- clone.files.in.dir[order(as.numeric(gsub(".*_S|_.*$|\\..*$", '', clone.files.in.dir)))]
count.files.in.dir <- list.files(count.dir);
count.files.in.dir <- count.files.in.dir[order(as.numeric(gsub(".*_S|\\..*", '', count.files.in.dir)))]

if (length(clone.files.in.dir) != length(count.files.in.dir)) stop("Files do not match")


# Create empty arrays
raw.entropies <- numeric(length(clone.files.in.dir)); naive.entropies <- raw.entropies; nb.entropies <- raw.entropies
unique.clones <- NULL
raw.clonality <- NULL; naive.clonality <- NULL; nb.clonality <- NULL
raw.max.clonal.freq <- numeric(length(count.files.in.dir)); naive.max.clonal.freq <- raw.max.clonal.freq; nb.max.clonal.freq <- raw.max.clonal.freq
raw.norm.entropy <- numeric(length(clone.files.in.dir)); naive.norm.entropy <- raw.norm.entropy; nb.norm.entropy <- raw.norm.entropy
raw.adaptive.clonality <- NULL; naive.adaptive.clonality <- raw.adaptive.clonality; nb.adaptive.clonality <- raw.adaptive.clonality
raw.max.clone.count <- numeric(length(clone.files.in.dir)); naive.max.clone.count <- raw.max.clone.count; nb.max.clone.count <- raw.max.clone.count
raw.gini <- numeric(length(clone.files.in.dir)); nb.gini <- numeric(length(clone.files.in.dir)); naive.gini <- nb.gini
raw.true <- numeric(length(clone.files.in.dir)); nb.true <- numeric(length(clone.files.in.dir)); naive.true <- nb.true

### Set divisions, if used
if (!(is.na(divisions_v)) & ((divisions_v == TRUE) | (divisions_v == T))) {
    divisions_lslsmat <- list()
    types_v <- c("raw", "naive", "nb")
    toDivide_v <- c(10, 25, 50, 100, 200, 250, 500)
    for (j in 1:length(types_v)){
        currType_v <- types_v[j]
        divisions_lslsmat[[currType_v]] <- list()
        for (i in 1:length(toDivide_v)){
            currDiv_v <- toDivide_v[i]
            currMat_mat <- data.frame(matrix(nrow = currDiv_v, ncol = length(clone.files.in.dir)))
            currName_v <- paste(currType_v, "top", currDiv_v, sep = '.')
            divisions_lslsmat[[currType_v]][[currName_v]] <- currMat_mat
        } # for i
    } # for
} # fi

        
    


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

    ## Create names variable, in case it's used later
    names_v <- c("Clone fraction", column, "nb.clone.fraction")

    ###
    ### CALCULATIONS
    ###

    ## count number of lines in file, i.e. number of unique clonotypes
    unique.clones[i] <- clone.curr.record[,.N]

    ##   calculate entropy
    raw.entropies[i] <- entropy::entropy(clone.curr.record[["Clone fraction"]], method="ML", unit="log");
    naive.entropies[i] <- entropy::entropy(clone.curr.record[[column]], method="ML", unit="log");
    nb.entropies[i] <- entropy::entropy(clone.curr.record[["nb.clone.fraction"]], method="ML", unit="log");
    

    ##   calculate clonality
    raw.clonality[i] <- 1 - (raw.entropies[i] / log(unique.clones[i]))
    naive.clonality[i] <- 1 - (naive.entropies[i] / log(unique.clones[i]))
    nb.clonality[i] <- 1 - (nb.entropies[i] / log(unique.clones[i]))

    ##   calculate clonality as the inverse of normalized shannon entropy
    ## Normalized shannon entropy
    
    raw.norm.entropy[i] <- raw.entropies[i] / log(unique.clones[i])
    naive.norm.entropy[i] <- naive.entropies[i] / log(unique.clones[i])
    nb.norm.entropy[i] <- nb.entropies[i] / log(unique.clones[i])

    ## New clonality
    raw.adaptive.clonality[i] <- 1 / raw.norm.entropy[i]
    naive.adaptive.clonality[i] <- 1 / naive.norm.entropy[i]
    nb.adaptive.clonality[i] <- 1 / nb.norm.entropy[i]

    ## tcR Gini Index
    raw.gini <- gini(.data = clone.curr.record[["Clone fraction"]], .do.norm = F)
    nb.gini <- gini(.data = clone.curr.record[["nb.clone.fraction"]], .do.norm = F)
    naive.gini <- gini(.data = clone.curr.record[[column]], .do.norm = F)

    ## tcR True Diversity
    raw.true <- diversity(.data = clone.curr.record[["Clone fraction"]], .do.norm = F)
    nb.true <- diversity(.data = clone.curr.record[["nb.clone.fraction"]], .do.norm = F)
    naive.true <- diversity(.data = clone.curr.record[[column]], .do.norm = F)
	

    ## Change clone frequency column to a percentage
    clone.curr.record[["Clone fraction"]] <- clone.curr.record[["Clone fraction"]] * 100
    clone.curr.record[[column]] <- clone.curr.record[[column]] * 100
    clone.curr.record[["nb.clone.fraction"]] <- clone.curr.record[["nb.clone.fraction"]] * 100
    
    
    ##	Calculate Max. clonotype frequency
    raw.max.clonal.freq[i] <- round(max(clone.curr.record[["Clone fraction"]]), digits = 4)
    naive.max.clonal.freq[i] <- round(max(clone.curr.record[[column]]), digits = 4)
    nb.max.clonal.freq[i] <- round(max(clone.curr.record[["nb.clone.fraction"]]), digits = 4)
    

    ##   Record maximum clone count
    raw.max.clone.count[i] <- max(clone.curr.record$`Clone count`)
    naive.max.clone.count[i] <- max(clone.curr.record$`Normalized clone count`)
    nb.max.clone.count[i] <- max(clone.curr.record$`nb.clone.count`)

    ## Record top clone values for each normalization method
    if (!(is.na(divisions_v)) & ((divisions_v == TRUE) | (divisions_v == T))) {
        for (j in 1:length(types_v)){
            ## Get norm type and associated mixcr column name
            currType_v <- types_v[j]
            currName_v <- names_v[j]
            ## Sort clones by that column
            clone.curr.record <- clone.curr.record[order(clone.curr.record[[currName_v]], decreasing = T),]
            ## Get values for each top clone divisions
            for (k in 1:length(toDivide_v)){
                currDivide_v <- toDivide_v[k]
                currDivName_v <- paste(currType_v, "top", currDivide_v, sep = '.')
                divisions_lslsmat[[ currType_v ]][[ currDivName_v ]][,i] <- clone.curr.record[[ currName_v ]][1:currDivide_v]
            } # for k
        } # for j
    } # fi
    


    #   update progress
    if((i %%10) == 0)   {
        cat("Processing file ", i, " (out of ", length(count.files.in.dir), ")\n", sep="");
    }   #   fi

}	#	for i

if (!(is.na(divisions_v)) & ((divisions_v == TRUE) | (divisions_v == T))) {
    for (i in 1:length(divisions_lslsmat)){
        currNorm_lsmat <- divisions_lslsmat[[i]]
        for (j in 1:length(currNorm_lsmat)){
            currDiv_mat <- currNorm_lsmat[[j]]
            currSummary_mat <- round(t(apply(currDiv_mat, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
            colnames(currSummary_mat) <- paste(names(currNorm_lsmat)[j], c("mean", "median", "sum"), sep = '.')
            currNorm_lsmat[[j]] <- currSummary_mat
        } # for j
        divisions_lslsmat[[i]] <- currNorm_lsmat
    } # for i
} # fi



#   create output data.frame
output.df <- data.frame(clone.files.in.dir, unique.clones, raw.entropies, naive.entropies, nb.entropies,
                        raw.norm.entropy, naive.norm.entropy, nb.norm.entropy,
                        raw.clonality, naive.clonality, nb.clonality,
                        raw.adaptive.clonality, naive.adaptive.clonality, nb.adaptive.clonality,
	     		raw.max.clonal.freq, naive.max.clonal.freq, nb.max.clonal.freq,
                        raw.max.clone.count, naive.max.clone.count, nb.max.clone.count,
                        raw.gini, naive.gini, nb.gini, raw.true, naive.true, nb.true)


outCols_v <- c("File", "Unique Clonotypes", "raw.shannon", "naive.shannon", "nb.shannon",
               "raw.norm.shannon", "naive.norm.shannon", "nb.norm.shannon",
               "raw.Clonality", "naive.Clonality", "nb.Clonality",
               "raw.Adaptive.Clonality", "naive.Adaptive.Clonality", "nb.Adaptive.Clonality",
               "raw.Max.Clonal.Freq", "naive.Max.Clonal.Freq", "nb.Max.Clonal.Freq",
               "raw.Max.Clone.Count", "naive.Max.Clone.Count", "nb.Max.Clone.Count",
               "raw_gini", "naive_gini", "nb_gini", "raw_trueD", "naive_trueD", "nb_trueD")

colnames(output.df) <- outCols_v

if (!(is.na(divisions_v)) & ((divisions_v == TRUE) | (divisions_v == T))) {

    raw_df <- do.call("cbind", divisions_lslsmat$raw)
    naive_df <- do.call("cbind", divisions_lslsmat$naive)
    nb_df <- do.call("cbind", divisions_lslsmat$nb)

    output.df <- data.frame(output.df, raw_df, naive_df, nb_df)
} # fi


#   write output
file.name <- "tcR.multi.uniques.shannon.clonality.txt"
write.table(output.df, 
            file=paste(out.dir, file.name, sep = ''),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
