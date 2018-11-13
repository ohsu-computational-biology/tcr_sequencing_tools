#!/usr/bin/Rscript

### Script to perform a homeostatic analysis based off of the anlaysis in the tcR package.
### For a given treatment group, take all of the samples and for each sample:
### Find the average frequency of each clone in a certain frequency group (Rare, Small, medium, etc.)
### Also find the total frequency of all clones in that group (this is what the homeostatic analysis does)

### Output results for both mean frequency as well as cumulative frequency

    ### Mean Freq - take the average frequency of all of the clones in a particular division
	### For example, a sample may have 10k "Rare" clones, all with frequencies between 0 and 0.00001
	### This means their average frequency will also be somewhere within that range and gives an idea
	### of what the "typical" Rare clone looks like. These results are not often used.

    ### Cumulative Freq - take the sum of all of the frequencies of all of the clones in a particular division
	### For example, that same sample of 10k "Rare" clones, when all of their frequencies are summed,
	### will almost certainly have a value greater than 0.00001. The 'cumulative freq' is the 'homeostatic' result.
	### This result tells us what proportion of all of the clones are considered Rare clones.

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
suppressMessages(library(optparse))
source("/home/exacloud/lustre1/CompBio/users/hortowe/2016_11_27_stable_repos/WesPersonal/utilityFxns.R")

### Command LIne
optlist <- list(
	make_option(
		c("-c", "--cloneDir"),
		type = "character",
		help = "Path to directory of normalized clone files."
	),
	make_option(
		c("-m", "--metadata"),
		type = "character",
		help = "Path to metadata file."
	),
	make_option(
		c("-o", "--outDir"),
		type = "character",
		help = "Path to output directory"
	),
	make_option(
		c("-w", "--write"),
		type = "logical",
		help = "TRUE - write output. FALSE - no output"
	),
	make_option(
		c("-l", "--old"),
		type = "logical",
		help = "TRUE - old column names; FALSE - new column names"
	),
	make_option(
		c("-t", "--tissue"),
		type = "character",
		help = "Specific tissue to subset by. If blank, will use all tissues"
	),
	make_option(
		c("-y", "--type"),
		type = "character",
		help = "Character vector for a specific category of treatments to divide by, rather than each treatment individually.
			If blank, will not divide."
	)
)

### Parse Command Line
p <- OptionParser(usage = "%prog -c cloneDir -m metadata -o outDir -w write -l old -t tissue -y type",
		option_list = optlist)
args <- parse_args(p)
opt <- args$options

### Get command-line arguments
cloneDir_v <- args$cloneDir
metadata_v <- args$metadata
outDir_v <- args$outDir
toWrite_v <- args$write
old_v <- args$old
tissue_v <- args$tissue
type_v <- args$type

print("After assign args")
### Get files and names
cloneFiles_v <- list.files(cloneDir_v)
cloneFiles_v <- cloneFiles_v[order(as.numeric(gsub("^.*_S|_.*$", "", cloneFiles_v)))]
cloneNames_v <- sapply(cloneFiles_v, function(x) {
        temp <- unlist(strsplit(x, split = "_"))
        name_v <- grep("S[0-9]+", temp, value = T)
        return(name_v)}, USE.NAMES = F)
names(cloneFiles_v) <- cloneNames_v
batchName_v <- strsplit(cloneFiles_v[1], split = "_")[[1]][1]

### Read in metadata
metadata_dt <- fread(metadata_v)

### Extract samples from metadata
sampleCol_v <- grep("[Ss]ample", colnames(metadata_dt), value = T)[1]
samples_v <- paste0("S", gsub("^S", "", metadata_dt[[sampleCol_v]]))

### Subset cloneFiles_v and cloneNames_v for these samples
cloneNames_v <- intersect(samples_v, cloneNames_v)
cloneFiles_v <- cloneFiles_v[cloneNames_v]

### Read in metadata
metadata_dt <- fread(metadata_v)

### Read in data
cloneData_lsdt <- sapply(cloneFiles_v, function(x) fread(file.path(cloneDir_v, x)), simplify = F)

### Get fraction/count columns
## Get fraction column
if (old_v) {
    column_v <- "Normalized clone fraction"
    count_v <- "Normalized clone count"
} else {
    column_v <- "nb.clone.fraction"
    count_v <- "nb.clone.count"
} # fi

### Change columns in case it's raw data
if (!column_v %in% colnames(cloneData_lsdt[[1]])){
    column_v <- grep("cloneFraction|Clone fraction", colnames(cloneData_lsdt[[1]]), value = T)
    count_v <- grep("cloneCount|Clone count", colnames(cloneData_lsdt[[1]]), value = T)
} # fi

### Get appropriate columns and subset
tissueCol_v <- grep("issue", colnames(metadata_dt), value = T)

if (!is.null(tissue_v)) {
    ## Subset metadata
    metadata_dt <- metadata_dt[get(tissueCol_v) == tissue_v,]
    ## Subset clone data
    keepSamp_v <- paste0("S", metadata_dt[, get(sampleCol_v)])
    cloneData_lsdt <- cloneData_lsdt[keepSamp_v]
} #fi

### Get treatments to run analysis on
if (!is.null(type_v)){
    treatCol_v <- type_v
} else {
    treatCol_v <- grep("[Tt]reatment", colnames(metadata_dt), value = T)
} # fi

treatments_v <- unique(metadata_dt[, get(treatCol_v)])

### Create divisions and empty data.tables to hold results.
divisions_v <- c("Blank" = 0, "Rare" = 0.00001, "Small" = 0.0001, "Medium" = 0.001, "Large" = 0.01, "Hyperexpanded" = 1)

### Per-Sample Summary tables
    ### One entry for each table recording meanFreq, cumFreq, or # clones per freqGroup
    ### Simply concatenating the individual treatment tables together, for easier comparison.
perSampleMeanFreq_dt <- NULL
perSampleCumFreq_dt <- NULL
perSampleGrpCount_dt <- NULL

### Per-Treatment Summary tables
    ### One entry for each treatment. Takes the mean of the per-sample results.
    ### end up with mean(meanFreq); mean(cumFreq), mean(# clones/group)
perTreatMeanFreq_dt <- NULL
perTreatCumFreq_dt <- NULL
perTreatGrpCount_dt <- NULL

if (toWrite_v){
    ## file for graphs
    pdf(file = paste0(outDir_v, batchName_v, "_cumFreqHomeo.pdf"))
    ## Directories
    meanDir_v <- mkdir(outDir_v, "meanFreq")
    cumDir_v <- mkdir(outDir_v, "cumFreq")
} # fi

for (i in 1:length(treatments_v)){

    ## Get Treatment
    currTreat_v <- treatments_v[i]

    ## Subset files
    currFiles_v <- paste0("S", metadata_dt[get(treatCol_v) == currTreat_v, get(sampleCol_v)]) # LIB170728LC
    currData_lsdt <- cloneData_lsdt[currFiles_v]

    ## Make empty matrix
    meanFreq_mat <- matrix(ncol = (length(divisions_v)-1), nrow = length(currData_lsdt))
    cumFreq_mat <- meanFreq_mat
    groupCount_mat <- meanFreq_mat     # Record number of clones in each freq group for each treatment

    ## Print nrow (n clones) of each sample
    nrow_v <- sapply(currData_lsdt, nrow)
    print(c("Treatment" = currTreat_v, nrow_v))

    ## Remove zero-count clones
    currData_lsdt <- lapply(currData_lsdt, function(x) x[get(column_v) > 0,])

    ## Compare
    newNrow_v <- sapply(currData_lsdt, nrow)
    diff_v <- nrow_v - newNrow_v
    print("Removed the following empty clones from each sample:")
    print(c("Treatment" = currTreat_v, diff_v))

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
            ## get mean Freq and cumulative freq
            currMeanFreq_v <- mean(currSub_dt[,get(column_v)])
            currCumFreq_v <- sum(currSub_dt[,get(column_v)])
            ## Update matrices
            meanFreq_mat[k,(j-1)] <- currMeanFreq_v
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
        ## Write mean freq
        write.table(meanFreq_dt, file = file.path(meanDir_v, paste0(currTreat_v, "_meanClonalFreqs.txt")),
                    sep = '\t', quote = F, row.names = F)
        ## Write cum freq
        write.table(cumFreq_dt, file = file.path(cumDir_v, paste0(currTreat_v, "_cumClonalFreqs.txt")),
                    sep = '\t', quote = F, row.names = F)
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

    ## Plot
    finalHomeo_gg <- vis.clonal.space(perTreatCumFreq_mat) +
        ggtitle(paste0("Treat-wise Clonal Space Homeostasis (Cum. Freq.)")) +
        theme(plot.title = element_text(hjust = 0.5))
    print(finalHomeo_gg)
    dev.off()
} # fi
