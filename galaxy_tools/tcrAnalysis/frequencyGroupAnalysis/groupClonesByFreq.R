#!/usr/bin/Rscript

###
### Group Clones by different frequency classifications
###

### Take a batch of clone files, combine all of the files togther, and group them by their different frequencies

####################
### DEPENDENCIES ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

suppressMessages(library(data.table))

#################
### FUNCTIONS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################

mergeDTs <- function(data_lsdt, mergeCol_v, keepCol_v = NULL, ...) {
  #' Merge many data.tables together
  #' @description Take many data.tables and merge on and ID column, extracting a single column from each data.table as the column of interest
  #' @param data_lsdt list of data.tables to merge
  #' @param mergeCol_v which column from all of the data.tables to use to merge
  #' @param keepCol_v which column from all of the data.tables to use as the column of interest. If NULL, use all columns
  #' @param ... extra parameters passed to merge
  #' @value data.table with ncol == length(data_lsdt) + 1. Column names are names of list, or defaults to V1, V2,...
  #' @export
  
  ## Grab extra arguments
  extraParams_lsv <- list(...)
  
  ## Handle extra arguments
  if (!is.null(extraParams_lsv$all)){
    all_v <- extraParams_lsv$all
  } else {
    all_v <- T
  } # fi
  
  if (!is.null(extraParams_lsv$sort)){
    sort_v <- extraParams_lsv$sort
  } else {
    sort_v <- F
  } # fi
  
  ## If keepCol_v is NULL, grab all other columns
  if (is.null(keepCol_v)){
    keepCol_v <- colnames(data_lsdt[[1]])[-which(colnames(data_lsdt[[1]]) %in% mergeCol_v)]
  } # fi
  
  ## Create initial table by extracting the 2 columns of interest from the rest
  merge_dt <- data_lsdt[[1]][,mget(c(mergeCol_v, keepCol_v))]
  
  ## Create initial column names (first check if list has names and add if not)
  if (is.null(names(data_lsdt))) {
    names_v <- paste("V", 1:length(data_lsdt))
    names(data_lsdt) <- names_v
  } # fi
  
  if (length(keepCol_v) > 1){
    colNames_v <- c(mergeCol_v, paste(names(data_lsdt)[1], keepCol_v, sep = "_"))
  } else {
    colNames_v <- c(mergeCol_v, names(data_lsdt)[1])
  } # fi
  
  for (i in 2:length(data_lsdt)) {
    merge_dt <- merge(merge_dt,
                      data_lsdt[[i]][,mget(c(mergeCol_v, keepCol_v))],
                      by = mergeCol_v,
                      all = all_v, sort = sort_v)
    ## Update column names
    if (length(keepCol_v) > 1){
      colNames_v <- c(colNames_v, paste(names(data_lsdt)[i], keepCol_v, sep = "_"))
    } else {
      colNames_v <- c(colNames_v, names(data_lsdt)[i])
    } # fi
    
    ## Rename columns
    colnames(merge_dt) <- colNames_v
  } # for
  return(merge_dt)
}

####################
### COMMAND LINE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

### Get command-line arguments
arguments <- commandArgs(trailingOnly = T)

inputFiles_v <- arguments[1]
inputNames_v <- arguments[2]
columns_v <- arguments[3] ### Columns to read in ("nb.clone.count,nb.clone.fraction,clonalSequence,aaSeqCDR3,V segments,J segments")
meta_v <- arguments[4]
rare_v <- arguments[5]
small_v <- arguments[6]
medium_v <- arguments[7]
large_v <- arguments[8]
hyper_v <- arguments[9]
full_v <- arguments[10]
summary_v <- arguments[11]
outputs_lsv <- list("Rare" = rare_v, "Small" = small_v, "Medium" = medium_v, "Large" = large_v, 
                    "Hyperexpanded" = hyper_v, "full" = full_v, "summary" = summary_v)

### Print
print(inputFiles_v)
print(inputNames_v)
print(columns_v)
print(meta_v)

### Handle non-traditional arguments
inputFiles_v <- unlist(strsplit(inputFiles_v, split = ","))
columns_v <- unlist(strsplit(columns_v, split = ','))

### Divisions
divisions_v <- c(0.00001, 0.0001, 0.001, 0.01, 1)
names(divisions_v) <- c("Rare", "Small", "Medium", "Large", "Hyperexpanded")

#############
### INPUT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############

### Read in names
inputNames_dt <- fread(inputNames_v, header = F)

### Read metadata
meta_dt <- fread(meta_v)

### Get sample, treatment, and clone freq columns
sampleCol_v <- grep("ample", colnames(meta_dt), value = T)
treatCol_v <- grep("eatment", colnames(meta_dt), value = T)
column_v <- grep("fraction", columns_v, value = T)
if (length(column_v) == 0) column_v <- "normFreq" # If batch-normalized data, frequency column will be different

### NOT SURE HOW THIS WORKS BECAUSE batch_v will likely have length > 1, but grep doesn't work like that...
### Pre-subset inputFiles for specific batch, if necessary
# batchCol_v <- grep("atch", colnames(meta_dt), value = T)
# if (length(batchCol_v) == 1) {
#     batch_v <- unique(meta_dt[[batchCol_v]])
#     inputFiles_v <- grep(batch_v, inputFiles_v, value = T)
#     print(sprintf("Pre-subset inputFiles_v to only contain files from batch %s", batch_v))
# }

###
### Subset to only contain samples in meta
### Only works if sample column is # or S#
###

### Get sample numbers
samples_v <- gsub("^S", "", meta_dt[[sampleCol_v]])

### Get prefix and suffix of input names, then combine with sample numbers
baseFile_v <- strsplit(inputNames_dt$V1[1], split = "S[0-9]+")[[1]]
toKeep_v <- paste0(baseFile_v[1], "S", samples_v, baseFile_v[2])

### Get indeces of where they are in the inputNames data.table
inputIDX_v <- which(inputNames_dt$V1 %in% toKeep_v)

### Subset each
inputFiles_v <- inputFiles_v[inputIDX_v]
inputNames_dt <- inputNames_dt[inputIDX_v,]

### Read in data
clones_lsdt <- NULL
for (i in 1:length(inputFiles_v)) {
  ## Get data
  currData_dt <- fread(inputFiles_v[i], select = columns_v)
  currName_v <- inputNames_dt$V1[i]
  ## Get sample from input name
  currSample_v <- grep("S[0-9]+", unlist(strsplit(currName_v, split = "_")), value = T)
  currSampNum_v <- as.numeric(gsub("S", "", currSample_v))
  ## Get treatment from metadata
  currTreat_v <- meta_dt[get(sampleCol_v) %in% c(currSample_v, currSampNum_v), get(treatCol_v)]
  ## Add Sample and treatment to data. Also add clone index column and empty column for divisions
  currData_dt$Sample <- currSample_v
  currData_dt$Treatment <- currTreat_v
  currData_dt$id <- 1:nrow(currData_dt)
  currData_dt$Div <- character()
  ## Add to list
  clones_lsdt[[currName_v]] <- currData_dt
}

print("Finished reading in data.")

### Remove any empty files
clones_lsdt <- clones_lsdt[sapply(clones_lsdt, function(x) dim(x)[1]) > 0]

### Remove any empty clones
clones_lsdt <- lapply(clones_lsdt, function(x) x[get(column_v) > 0,])

### Classify Clones
clones_lsdt <- sapply(clones_lsdt, function(x) {
    ## Iterate for each row and determine which divisions are greater than the current row's frequency
    ## Since the division values represent the upper limit of that grouping, the first division that is
    ## greater than the frequency is the one that we want. (e.g. Large is 0.001 to 0.01, Medium is 0.0001 to 0.001,
    ## and Hyper is 0.01 to 1. Given freq of 0.00099, which is less than all 3, medium will be chosen. Given freq of
    ## freq of 0.0010001, which is less than Large and Hyper, Large will be chosen. Given freq of 0.01, Large and Hyper
    ## will match, and large will be chosen.)
    x[, Div := sapply(x[,get(column_v)], function(y) names(which(divisions_v >= y)[1]))]
    return(x)
}, simplify = F)

### Combine into giant data.table
clones_dt <- do.call(rbind, clones_lsdt)
print("Finished classification.")

### Gather summary info
summary_lsdt <- sapply(names(divisions_v), function(x) clones_dt[Div == x, .N, by = Sample], simplify = F)
summary_dt <- mergeDTs(summary_lsdt, mergeCol_v = "Sample", keepCol_v = "N")

### Split into multiple data.tables
divisions_lsdt <- list("summary" = summary_dt, "full" = clones_dt)
for (i in 1:length(divisions_v)){
    ## Get divisions
    currDiv_v <- names(divisions_v)[i]
    ## Subset
    currSubset_dt <- clones_dt[Div == currDiv_v,]
    ## Add to list
    divisions_lsdt[[currDiv_v]] <- currSubset_dt
} # for i
print("Finish split and summary.")

### Write out all of the data.tables
for (i in 1:length(divisions_lsdt)){
  ## Get division  name
  currName_v <- names(divisions_lsdt)[i]
  
  ## Get data
  currOut_dt <- divisions_lsdt[[currName_v]]
    
  ## Get output file
  currOut_v <- outputs_lsv[[currName_v]]
  
  ## Write table
  write.table(currOut_dt,
              file = currOut_v,
              sep = '\t', quote = F, row.names = F)
} # for i
print("Finish writing.")

warnings()

