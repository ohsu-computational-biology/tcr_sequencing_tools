#!/usr/bin/Rscript

### Script to perform a homeostatic analysis based off of the anlaysis in the tcR package.
### For a given treatment group, take all of the samples and for each sample:
### Find the average frequency of each clone in a certain frequency group (Rare, Small, medium, etc.)
### Also find the total frequency of all clones in that group (this is what the homeostatic analysis does)

### Requires
	### directory to normalized clone count files
	### metadata file with sample names and treatments
	### path to output
	### tissue to subset analysis by (optional)
	### treatment to subset by (optional)
	### TRUE/FALSE if pdf and text output should be written
### Note, this does everything that the homeoAnalysis.R does, plus more.

###   Load necessary libraries
print("Start")
suppressMessages(library(entropy))
suppressMessages(library(data.table))
suppressMessages(library(tcR))
suppressMessages(library(ggplot2))
suppressMessages(library(xlsx))
source("/home/exacloud/lustre1/CompBio/users/hortowe/2016_11_27_stable_repos/WesPersonal/utilityFxns.R")

###	Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);

cloneDir_v <- arguments[1];    # Typically .../dhaarini/DNAXXXXLC/normalization/normalized_clones/
metadata_v <- arguments[2]     # Typically $data/QC/meta
outDir_v <- arguments[3];      # Typeically $data/QC/clonalFreq
toWrite_v <- arguments[4]
tissue_v <- arguments[5]       # specify a specific tissue to subset by. If NA, will use all tissues
type_v <- arguments[6]         # sometimes divide by a certain category of treatment, rather than each treatment. If NA, will not divide

print("After assign args")
### Get files
cloneFiles_v <- list.files(cloneDir_v)
cloneFiles_v <- cloneFiles_v[order(as.numeric(gsub("^.*_S|_.*$", "", cloneFiles_v)))]

### Get batches and prefix
batches_v <- unique(gsub("_S.*$", "", cloneFiles_v))
pre_v <- gsub("[0-9]*LC", "", batches_v[1])

### Get sample names in format of batch_S#
cloneNames_v <- sapply(cloneFiles_v, function(x) {
        temp <- unlist(strsplit(x, split = "_"))
        sample_v <- grep("S[0-9]+", temp, value = T)
	batch_v <- temp[1]
	name_v <- paste(batch_v, sample_v, sep = "_")
        return(name_v)}, USE.NAMES = F)

### Read in data
cloneData_lsdt <- sapply(cloneFiles_v, function(x) fread(file.path(cloneDir_v, x)), simplify = F)
names(cloneData_lsdt) <- cloneNames_v

metadata_dt <- fread(metadata_v)

### Get fraction/count columns
column_v <- "normFreq"
count_v <- "normC"

### Get appropriate columns and subset

tissueCol_v <- grep("issue", colnames(metadata_dt), value = T)
sampleCol_v <- grep("ample", colnames(metadata_dt), value = T)[1]
batchCol_v <- grep("atch", colnames(metadata_dt), value = T)

if (!is.na(tissue_v)) {
    ## Subset metadata
    metadata_dt <- metadata_dt[get(tissueCol_v) == tissue_v,]
    ## Subset clone data
    keepSamp_v <- paste0("S", metadata_dt[, get(sampleCol_v)])
    ## Add batch to sample name, if needed
    if (!keepSamp_v[1] %in% names(cloneData_lsdt)) {
	keepSamp_v <- paste0(pre_v, metadata_dt[,get(batchCol_v)], "LC_S", metadata_dt[,get(sampleCol_v)])
    } # fi
    ## Subset
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

perSampleMeanFreq_dt <- NULL
perSampleCumFreq_dt <- NULL
perSampleGrpCount_dt <- NULL
perTreatMeanFreq_dt <- NULL
perTreatCumFreq_dt <- NULL
perTreatGrpCount_dt <- NULL

if (toWrite_v){
    ## file for graphs
    if (!is.na(tissue_v)) batchName_v <- tissue_v

    pdf(file = file.path(outDir_v, paste0(batchName_v, "_cumFreqHomeo.pdf")))
    ## Directories
    meanDir_v <- mkdir(outDir_v, "meanFreq")
    cumDir_v <- mkdir(outDir_v, "cumFreq")
} # fi

for (i in 1:length(treatments_v)){

    ## Get Treatment
    currTreat_v <- treatments_v[i]

    ## Subset files
    currFiles_v <- paste0("S", metadata_dt[get(treatCol_v) == currTreat_v, get(sampleCol_v)]) # LIB170728LC

    ## Add batch to sample name, if needed
    if (!currFiles_v[1] %in% names(cloneData_lsdt)) {
        currFiles_v <- paste0(pre_v, 
				metadata_dt[get(treatCol_v) == currTreat_v, get(batchCol_v)], 
				"LC_S", 
				metadata_dt[get(treatCol_v) == currTreat_v, get(sampleCol_v)])
    } # fi
    ## Subse
    currData_lsdt <- cloneData_lsdt[currFiles_v]

    ## Make empty matrix
    meanFreq_mat <- matrix(ncol = (length(divisions_v)-1), nrow = length(currData_lsdt))
    cumFreq_mat <- meanFreq_mat
    groupCount_mat <- meanFreq_mat

    ## Print nrow
    nrow_v <- sapply(currData_lsdt, nrow)
    print(c("Treatment" = currTreat_v, nrow_v))

    ## Remove zero-count clones
    currData_lsdt <- lapply(currData_lsdt, function(x) x[get(column_v) > 0,])

    ## For each sample, get all of the clones in that
    for (j in 2:length(divisions_v)){
        ## Get cut-offs
        currMin_v <- divisions_v[j-1]
        currMax_v <- divisions_v[j]
        ## Subset for each sample
        for (k in 1:length(currData_lsdt)){
            ## Get sample
            currData_dt <- currData_lsdt[[k]]
            ## Subset
            currSub_dt <- currData_dt[(get(column_v) > currMin_v & get(column_v) <= currMax_v),]
            ## get mean Freq
            currMeanFreq_v <- mean(currSub_dt[,get(column_v)])
            ## Update matrix
            meanFreq_mat[k,(j-1)] <- currMeanFreq_v
            ## get cumulative freq
            currCumFreq_v <- sum(currSub_dt[,get(column_v)])
            ## Update matrix
            cumFreq_mat[k,(j-1)] <- currCumFreq_v
            ## Get number of clones and add to matrix
            groupCount_mat[k,(j-1)] <- currSub_dt[,.N]
        } # for k
    } # for j

    ## Convert to data table
    meanFreq_dt <- as.data.table(meanFreq_mat)
    cumFreq_dt <- as.data.table(cumFreq_mat)
    groupCount_dt <- as.data.table(groupCount_mat)

    ## Add division names as column names
    colnames(meanFreq_dt) <- names(divisions_v)[2:length(divisions_v)]
    colnames(cumFreq_dt) <- colnames(meanFreq_dt)
    colnames(groupCount_dt) <- colnames(meanFreq_dt)

    ## Add column of sample names
    samples_v <- names(currData_lsdt)
    meanFreq_dt$Sample <- samples_v
    cumFreq_dt$Sample <- samples_v
    groupCount_dt$Sample <- samples_v

    ## Reorder
    meanFreq_dt <- meanFreq_dt[,c(ncol(meanFreq_dt),1:(ncol(meanFreq_dt)-1)), with = F]
    cumFreq_dt <- cumFreq_dt[,c(ncol(cumFreq_dt),1:(ncol(cumFreq_dt)-1)), with = F]
    groupCount_dt <- groupCount_dt[,c(ncol(groupCount_dt),1:(ncol(groupCount_dt)-1)), with = F]

    ## Write out treatment table and make plot
    if (toWrite_v) {

	## Names
	if (!is.na(tissue_v)) {
		meanName_v <- paste0(tissue_v, "_", currTreat_v, "_meanClonalFreqs.txt")
		meanExcel_v <- paste0(tissue_v, "_meanClonalFreqs.xlsx")
		cumName_v <- paste0(tissue_v, "_", currTreat_v, "_cumClonalFreqs.txt")
		cumExcel_v <- paste0(tissue_v, "_cumClonalFreqs.xlsx")
	} else {
		meanName_v <- paste0(currtreat_v, "_meanClonalFreqs.txt")
		meanExcel_v <- "meanClonalFreqs.xlsx"
		cumName_v <- paste0(currtreat_v, "_cumClonalFreqs.txt")
		cumExcel_v <- "cumClonalFreqs.xlsx"
	} # fi

        ## Write mean freq
        write.table(meanFreq_dt, file = file.path(meanDir_v, meanName_v),
                    sep = '\t', quote = F, row.names = F)
        ## Write cum freq
        write.table(cumFreq_dt, file = file.path(cumDir_v, cumName_v), 
                    sep = '\t', quote = F, row.names = F)

	## Write excel
	if (i == 1) {append_v <- F} else {append_v <- T}
	write.xlsx(meanFreq_dt, file = file.path(meanDir_v, meanExcel_v), sheetName = currTreat_v,
			row.names=F, append = append_v)
	write.xlsx(cumFreq_dt, file = file.path(cumDir_v, cumExcel_v), sheetName = currTreat_v,
			row.names=F, append = append_v)

        ## Make plot
        colnames(cumFreq_mat) <- names(divisions_v)[2:length(divisions_v)]
        rownames(cumFreq_mat) <- samples_v
        currHomeo_gg <- vis.clonal.space(cumFreq_mat) +
            ggtitle(paste0(currTreat_v, " Clonal Space Homeostasis (Cum. Freq.)")) +
            theme(plot.title = element_text(hjust = 0.5))
        print(currHomeo_gg)
    } # fi

    ## Combine with overall table
    perSampleMeanFreq_dt <- rbind(perSampleMeanFreq_dt, meanFreq_dt)
    perSampleCumFreq_dt <- rbind(perSampleCumFreq_dt, cumFreq_dt)
    perSampleGrpCount_dt <- rbind(perSampleGrpCount_dt, groupCount_dt)

    ## Take mean for treat-wise table
    currMeanMean_v <- c(currTreat_v, apply(meanFreq_mat, 2, mean))
    currCumMean_v <- c(currTreat_v, apply(cumFreq_mat, 2, mean))
    currCountMean_v <- c(currTreat_v, apply(groupCount_mat, 2, mean))

    ## Construct treatwise
    perTreatMeanFreq_dt <- rbind(perTreatMeanFreq_dt, currMeanMean_v)
    perTreatCumFreq_dt <- rbind(perTreatCumFreq_dt, currCumMean_v)
    perTreatGrpCount_dt <- rbind(perTreatGrpCount_dt, currCountMean_v)
    
} # for i

### Finish construction
colnames(perTreatMeanFreq_dt) <- c("Treatment", names(divisions_v)[2:length(divisions_v)])
colnames(perTreatCumFreq_dt) <- colnames(perTreatMeanFreq_dt)
colnames(perTreatGrpCount_dt) <- colnames(perTreatMeanFreq_dt)

### Revert cumulative freq back to matrix for ggplot
if (nrow(perTreatCumFreq_dt) == 1) {
    perTreatCumFreq_mat <- t(as.matrix(perTreatCumFreq_dt[,2:ncol(perTreatCumFreq_dt)]))
} else {
    perTreatCumFreq_mat <- as.matrix(perTreatCumFreq_dt[,2:ncol(perTreatCumFreq_dt)])
} # fi
rownames(perTreatCumFreq_mat) <- perTreatCumFreq_dt[,1]
perTreatCumFreq_mat <- apply(perTreatCumFreq_mat, c(1,2), as.numeric)

### Write out
if (toWrite_v) {
    ## Mean Freq
    if (!is.na(tissue_v)){
	batchName_v <- tissue_v
    } # fi

    write.table(perSampleMeanFreq_dt, file = file.path(meanDir_v, paste0(batchName_v, "_meanClonalFreqs.txt")),
                sep = '\t', quote = F, row.names = F)
    write.table(perTreatMeanFreq_dt, file = file.path(meanDir_v, paste0(batchName_v, "_treatWiseMeanClonalFreqs.txt")),
                sep = '\t', quote = F, row.names = F)

    ## Cum Freq
    write.table(perSampleCumFreq_dt, file = file.path(cumDir_v, paste0(batchName_v, "_cumClonalFreqs.txt")),
                sep = '\t', quote = F, row.names = F)
    write.table(perTreatCumFreq_dt, file = file.path(cumDir_v, paste0(batchName_v, "_treatWiseCumClonalFreqs.txt")),
                sep = '\t', quote = F, row.names = F)

    ## Group Counts
    write.table(perSampleGrpCount_dt, file = file.path(outDir_v, paste0(batchName_v, "_groupCounts.txt")),
                sep = '\t', quote = F, row.names = F)
    write.table(perTreatGrpCount_dt, file = file.path(outDir_v, paste0(batchName_v, "_treatWiseGroupCounts.txt")),
                sep = '\t', quote = F, row.names = F)

    ## Write excel
    
    write.xlsx(perSampleMeanFreq_dt, file = file.path(meanDir_v, meanExcel_v), sheetName = "bySample",
		row.names = F, append = T)
    rownames(perTreatMeanFreq_dt) <- NULL
    write.xlsx(perTreatMeanFreq_dt, file = file.path(meanDir_v, meanExcel_v), sheetName = "byTreat",
		row.names = F, append = T)
    write.xlsx(perSampleCumFreq_dt, file = file.path(cumDir_v, cumExcel_v), sheetName = "bySample",
		row.names = F, append = T)
    rownames(perTreatCumFreq_dt) <- NULL
    write.xlsx(perTreatCumFreq_dt, file = file.path(cumDir_v, cumExcel_v), sheetName = "byTreat",
		row.names = F, append = T)

    ## Plot
    finalHomeo_gg <- vis.clonal.space(perTreatCumFreq_mat) +
        ggtitle(paste0("Treat-wise Clonal Space Homeostasis (Cum. Freq.)")) +
        theme(plot.title = element_text(hjust = 0.5))
    print(finalHomeo_gg)
    dev.off()
} # fi
