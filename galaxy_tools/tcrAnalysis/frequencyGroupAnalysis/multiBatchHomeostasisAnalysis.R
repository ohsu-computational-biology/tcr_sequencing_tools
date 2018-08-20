#!/usr/bin/Rscript

###################
### DESCRIPTION ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################

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

####################
### DEPENDENCIES ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

suppressMessages(library(data.table))
suppressMessages(library(tcR))
suppressMessages(library(ggplot2))
suppressMessages(library(xlsx))
suppressMessages(library(optparse))

####################
### COMMAND LINE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

### Command LIne
optlist <- list(
  make_option(
    c("-f", "--files"),
    type = "character",
    help = "Normalized clone files"
  ),
  make_option(
    c("-n", "--names"),
    type = "character",
    help = "Clone file names"
  ),
  make_option(
    c("-s", "--sampleID"),
    type = "numeric",
    help = "Index of sample number from sample name, when split by '_'."
  ),
  make_option(
    c("-m", "--metadata"),
    type = "character",
    help = "Path to metadata file."
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
  ),
  make_option(
    c("-p", "--plot"),
    type = "character",
    help = "Name of homeostasis output plot."
  ),
  make_option(
    c("-b", "--batchIndex"),
    type = "numeric",
    help = "Index of batch in sample name, if split by '_'."
  )
  )

### Parse Command Line
p <- OptionParser(usage = "%prog -f files -n names -s sampleID -m metadata -l old -t tissue -y type -p plot -b batchIndex",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

### Get command-line arguments
inputFiles_v <- args$files
inputNames_v <- args$names
sampleID_v <- args$sampleID
metadata_v <- args$metadata
old_v <- args$old
tissue_v <- args$tissue
type_v <- args$type
plot_v <- args$plot
batchIndex_v <- args$batchIndex

#############
### INPUT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############

### Get files and names
inputFiles_v <- unlist(strsplit(inputFiles_v, split = ','))
inputNames_dt <- fread(inputNames_v, header = F)

### Get batches and order them
batches_v <- unique(sapply(inputNames_dt$V1, function(x) strsplit(x, split = "_")[[1]][batchIndex_v]))
batches_v <- sort(batches_v)

### Get prefix? Not sure the point yet
pre_v <- gsub("[0-9]*LC", "", batches_v[1])

### Get sample names in format of batch_S#
inputNames_v <- sapply(inputNames_dt$V1, function(x) {
        temp <- unlist(strsplit(x, split = "_"))
        sample_v <- grep("S[0-9]+", temp, value = T)
	batch_v <- temp[1]
	name_v <- paste(batch_v, sample_v, sep = "_")
        return(name_v)}, USE.NAMES = F)

### Sort files by batch and then by sample within batch
names(inputFiles_v) <- inputNames_v
orderedNames_v <- unname(unlist(sapply(batches_v, function(x) {
  y <- inputNames_v[grep(x, inputNames_v)]
  z <- y[order(as.numeric(gsub(".*_S|_.*$|\\..*$}", '', y)))]
  return(z)}, simplify = F)))
inputFiles_v <- inputFiles_v[orderedNames_v]

### Read in data and metadata
cloneData_lsdt <- sapply(inputFiles_v, function(x) fread(x), simplify = F)
names(cloneData_lsdt) <- names(inputFiles_v)
metadata_dt <- fread(metadata_v)

### Get fraction/count columns and various metadata columns
column_v <- "normFreq"
count_v <- "normC"
tissueCol_v <- grep("[Tt]issue", colnames(metadata_dt), value = T)
sampleCol_v <- grep("[Ss]ample", colnames(metadata_dt), value = T)[1]
batchCol_v <- grep("[Bb]atch", colnames(metadata_dt), value = T)

if (tissue_v != "NA") {
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
if (type_v != "NA"){
    treatCol_v <- type_v
} else {
    treatCol_v <- grep("[Tt]reatment", colnames(metadata_dt), value = T)
} # fi

### Subset treatments, if needed
treatments_v <- unique(metadata_dt[, get(treatCol_v)])

### Create divisions
divisions_v <- c("Blank" = 0, "Rare" = 0.00001, "Small" = 0.0001, "Medium" = 0.001, "Large" = 0.01, "Hyperexpanded" = 1)

###################
### PREP OUTPUT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################

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

### List of individual treatment results
cumFreqOut_lsdt <- meanFreqOut_lsdt <- NULL

### List for plots
plotList_lsgg <- list()

################
### ANALYSIS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################ 

for (i in 1:length(treatments_v)){

    ## Get Treatment
    currTreat_v <- treatments_v[i]

    ## Subset files
    currFiles_v <- paste0("S", gsub("^S", "", metadata_dt[get(treatCol_v) == currTreat_v, get(sampleCol_v)])) # LIB170728LC
    
    ## Add batch to sample name, if needed (sometimes samples are listed as "S1" and other times "LIBXXXXLC_S1")
    if (!currFiles_v[1] %in% names(cloneData_lsdt)) {
        currFiles_v <- paste0(pre_v, 
				metadata_dt[get(treatCol_v) == currTreat_v, get(batchCol_v)], 
				"LC_S", 
				metadata_dt[get(treatCol_v) == currTreat_v, get(sampleCol_v)])
    } # fi
    
    ## Subset
    currData_lsdt <- cloneData_lsdt[currFiles_v]
    
    ## Make empty matrix
    meanFreq_mat <- matrix(ncol = (length(divisions_v)-1), nrow = length(currData_lsdt))
    cumFreq_mat <- meanFreq_mat
    groupCount_mat <- meanFreq_mat

    ## Print nrow
    nrow_v <- sapply(currData_lsdt, nrow)
    print("Treatment:"); print(currTreat_v)
    print("Nrow:"); print(nrow_v)

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
    
    ## Add to lists
    cumFreqOut_lsdt[[currTreat_v]] <- cumFreq_dt
    meanFreqOut_lsdt[[currTreat_v]] <- meanFreq_dt

    ## Make plot
    colnames(cumFreq_mat) <- names(divisions_v)[2:length(divisions_v)]
    rownames(cumFreq_mat) <- samples_v
    currHomeo_gg <- vis.clonal.space(cumFreq_mat) +
        ggtitle(paste0(currTreat_v, " Clonal Space Homeostasis (Cum. Freq.)")) +
        theme(plot.title = element_text(hjust = 0.5))
    plotList_lsgg[[i]] <- currHomeo_gg

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

######################
### WRANGLE OUTPUT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################

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

### Prepare data common to both excel sheets
preList_lsdt <- list("perSampleGroupCount" = perSampleGrpCount_dt, "treatWiseGroupCount" = perTreatGrpCount_dt)

### Prepare specific data
cumFreqSummary_lsdt <- list("perSampleFreq" = perSampleCumFreq_dt, "treatWiseFreq" = perTreatCumFreq_dt)
meanFreqSummary_lsdt <- list("perSampleFreq" = perSampleMeanFreq_dt, "treatWiseFreq" = perTreatMeanFreq_dt)

### Combine
cumFreqOut_lsdt <- c(preList_lsdt, cumFreqSummary_lsdt, cumFreqOut_lsdt)
meanFreqOut_lsdt <- c(preList_lsdt, meanFreqSummary_lsdt, meanFreqOut_lsdt)

##############
### OUTPUT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############

### Write
for (i in 1:length(cumFreqOut_lsdt)) {
  ## Get name and data
  currName_v <- names(cumFreqOut_lsdt)[i]
  currData_dt <- cumFreqOut_lsdt[[currName_v]]
  ## Determine if create or append
  append_v <- ifelse(i == 1, F, T)
  ## Write
  write.xlsx(currData_dt,
             file = "./tempCum.xlsx",
             sheetName = currName_v,
             row.names = F,
             append = append_v)
} # for i

for (i in 1:length(meanFreqOut_lsdt)) {
  ## Get name and data
  currName_v <- names(cumFreqOut_lsdt)[i]
  currData_dt <- cumFreqOut_lsdt[[currName_v]]
  ## Determine if create or append
  append_v <- ifelse(i == 1, F, T)
  ## Write
  write.xlsx(currData_dt,
             file = "./tempMean.xlsx",
             sheetName = currName_v,
             row.names = F,
             append = append_v)
} # for i 

## Final Plot
finalHomeo_gg <- vis.clonal.space(perTreatCumFreq_mat) +
  ggtitle(paste0("Treat-wise Clonal Space Homeostasis (Cum. Freq.)")) +
  theme(plot.title = element_text(hjust = 0.5))

plotList_lsgg[[(length(plotList_lsgg)+1)]] <- finalHomeo_gg

print("Plot list length")
print(length(plotList_lsgg))

## Output plot
pdf(plot_v)
for (i in 1:length(plotList_lsgg)) {
  print(plotList_lsgg[[i]])
}
graphics.off()


