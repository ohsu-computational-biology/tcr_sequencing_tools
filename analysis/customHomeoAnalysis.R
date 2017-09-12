#	RScript that calculates number of unique clonotypes, Shannon entropy, and clonality for multiple clonotype files
#
#   Input is a directory containing one or more files with MiXCR output format
#   
#   Load necessary libraries
library(entropy);
library(data.table);
library(tcR);
library(ggplot2);

#	Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);
#   Directory should contain all and only .txt files from exportClones, post-normalization with the following
#   format:
#     Clone count   Clone fraction    Clonal sequence(s)    AA. Seq. CDR3   Best V Hit    Best J Hit    V segments
#     J segments    Normalized clone count    Normalized clone fraction

cloneDir_v <- arguments[1];    # Typically .../dhaarini/DNAXXXXLC/normalization/normalized_clones/
metadata_v <- arguments[2]     # Typically $data/QC/meta
outDir_v <- arguments[3];      # Typeically $data/QC/clonalFreq
tissue_v <- arguments[4]       # specify a specific tissue to subset by. If NA, will use all tissues
type_v <- arguments[5]         # sometimes divide by a certain category of treatment, rather than each treatment. If NA, will not divide

### Get files and names
cloneFiles_v <- list.files(cloneDir_v)
cloneFiles_v <- cloneFiles_v[order(as.numeric(gsub("^.*_S|_.*$", "", cloneFiles_v)))]
cloneNames_v <- sapply(cloneFiles_v, function(x) strsplit(x, split = "_")[[1]][2], USE.NAMES=F)
batchName_v <- strsplit(cloneFiles_v[1], split = "_")[[1]][1]

### Read in data
cloneData_lsdt <- sapply(cloneFiles_v, function(x) fread(file.path(cloneDir_v, x)), simplify = F)
names(cloneData_lsdt) <- cloneNames_v

metadata_dt <- fread(metadata_v)

### Get appropriate columns and subset

tissueCol_v <- grep("issue", colnames(metadata_dt), value = T)
sampleCol_v <- grep("ample", colnames(metadata_dt), value = T)[1]

if (!is.na(tissue_v)) {
    ## Subset metadata
    metadata_dt <- metadata_dt[get(tissueCol_v) == tissue_v,]
    ## Subset clone data
    keepSamp_v <- paste0("S", metadata_dt[, get(sampleCol_v)])
    cloneData_lsdt <- cloneData_lsdt[keepSamp_v]
} #fi

### Get treatments to run analysis on
if (!is.na(type_v)){
    treatCol_v <- type_v
} else {
    treatCol_v <- grep("eatment", colnames(metadata_dt), value = T)
} # fi

treatments_v <- unique(metadata_dt[, get(treatCol_v)])


divisions_v <- c("Blank" = 0, "Rare" = 0.00001, "Small" = 0.0001, "Medium" = 0.001, "Large" = 0.01, "Hyperexpanded" = 1)
totalFreq_dt <- NULL
treatFreq_dt <- NULL

for (i in 1:length(treatments_v)){

    ## Get Treatment
    currTreat_v <- treatments_v[i]

    ## Subset files
    currFiles_v <- paste0("S", metadata_dt[get(treatCol_v) == currTreat_v, get(sampleCol_v)]) # LIB170728LC
    currData_lsdt <- cloneData_lsdt[currFiles_v]

    ## Make empty matrix
    meanFreq_mat <-	matrix(ncol = (length(divisions_v)-1), nrow = length(currData_lsdt))

    ## Print nrow
    nrow_v <- sapply(currData_lsdt, nrow)
    print(c("Treatment" = currTreat_v, nrow_v))

    ## Remove zero-count clones
    currData_lsdt <- lapply(currData_lsdt, function(x) x[`Normalized clone fraction` > 0,])

    ## For each sample, get all of the clones in that
    col_v <- "Normalized clone fraction"
    for (j in 2:length(divisions_v)){
        ## Get cut-offs
        currMin_v <- divisions_v[j-1]
        currMax_v <- divisions_v[j]
        ## Subset for each sample
        for (k in 1:length(currData_lsdt)){
            ## Get sample
            currData_dt <- currData_lsdt[[k]]
            ## Subset
            currSub_dt <- currData_dt[(get(col_v) > currMin_v & get(col_v) <= currMax_v),]
            ## get mean Freq
            currMeanFreq_v <- mean(currSub_dt[,get(col_v)])
            ## Update matrix
            meanFreq_mat[k,(j-1)] <- currMeanFreq_v
        } # for k
    } # for j

    ## Convert to data table
    meanFreq_dt <- as.data.table(meanFreq_mat)

    ## Add names
    colnames(meanFreq_dt) <- names(divisions_v)[2:length(divisions_v)]
    samples_v <- names(currData_lsdt)
    meanFreq_dt$Sample <- samples_v
    meanFreq_dt <- meanFreq_dt[,c(ncol(meanFreq_dt),1:(ncol(meanFreq_dt)-1)), with = F]

    ## Write out treatment table
    write.table(meanFreq_dt, file = file.path(outDir_v, paste0(currTreat_v, "_meanClonalFreqs.txt")),
                sep = '\t', quote = F, row.names = F)

    ## Combine with overall table
    totalFreq_dt <- rbind(totalFreq_dt, meanFreq_dt)

    ## Take mean for treat-wise table
    currMean_v <- c(currTreat_v, apply(meanFreq_mat, 2, mean))
    treatFreq_dt <- rbind(treatFreq_dt, currMean_v)


    
} # for i

### Finish construction
colnames(treatFreq_dt) <- c("Treatment", names(divisions_v)[2:length(divisions_v)])

### Write out
write.table(totalFreq_dt, file = file.path(outDir_v, paste0(batchName_v, "_meanClonalFreqs.txt")),
            sep = '\t', quote = F, row.names = F)
write.table(treatFreq_dt, file = file.path(outDir_v, paste0(batchName_v, "_treatWiseMeanClonalFreqs.txt")),
            sep = '\t', quote = F, row.names = F)

