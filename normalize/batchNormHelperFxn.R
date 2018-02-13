###
### Helper functions for batch normalization
###

### Read in a clone file
read_clone<-function(file, directory, norm_v) {
  #' Read in a clone file
  #' @description
  #' Designed to be used in conjunction with readDir_clone. Reads pertinent columns from a clone file into a data.table
  #' @param file Character string of the full name of the file to be read in
  #' @param directory path to the directory that contains 'file'
  #' @param norm_v character indicating if the raw clone count ("raw"), original norm counts ("orig"), or negative binomial norm count ("nb") should be used.
  #' @return data.table containing the Clone ID, Raw/Norm Clone Count, Clonal Sequence, V segment, and J segment information for an entire sample
  #' @export
  
  ## Specify desired columns
  if (norm_v == "orig"){
    getCols_v <- c("Clone ID", "cloneId", "Normalized clone count", "clonalSequence", "Clonal sequence(s)", "AA. Seq. CDR3", "aaSeqCDR3", "V segments", "J segments")
    newCols_v <- c("ID", "Seq", "aaSeq", "V", "J", "CloneCount")
  } else if (norm_v == "raw") {
    getCols_v <- c("Clone ID", "cloneId", "Clone count", "clonalSequence", "Clonal sequence(s)", "AA. Seq. CDR3", "aaSeqCDR3", "V segments", "J segments")
    newCols_v <- c("ID", "CloneCount", "Seq", "aaSeq", "V", "J")
  } else if (norm_v == "nb") {
    getCols_v <- c("Clone ID", "cloneId", "nb.clone.count", "clonalSequence", "Clonal sequence(s)", "AA. Seq. CDR3", "aaSeqCDR3", "V segments", "J segments")
    newCols_v <- c("ID", "Seq", "aaSeq", "V", "J", "CloneCount")
  } else {
    stop("Incorrect designation of norm_v argument. Must be 'orig', 'raw', or 'nb'.")
  }
  
  ## Read in first line (aka column names)
  colNames_v <- unlist(fread(file.path(directory, file), nrows = 1, sep = '\t', header = F))
  
  ## Get column indices (for some reason, fread doesn't work with character vector of column names, 
  ##                     so giving int vector of column indices)
  whichCols_v <- which(colNames_v %in% getCols_v)
  
  ## Handle edge case - Clone ID is repeated...need to remove second instance of it
  ## Making as general as possible, in case of different (or multiple) repeated columns
  colFreq_df <- as.data.frame(table(colNames_v[whichCols_v]))                # Get frequencies
  dupCols_v <- as.character(colFreq_df[colFreq_df$Freq > 1,"Var1"])          # Get repeated columns
  whichDup_v <- sapply(dupCols_v, function(x) {                              # Find whichCols_v indices that should be removed
    currCols_v <- which(colNames_v[whichCols_v] %in% x)
    badCols_v <- currCols_v[-1]
  })
  whichCols_v <- whichCols_v[-whichDup_v]
  
  ## Read in clone file
  clone<-fread(file.path(directory,file), select = whichCols_v, sep = '\t', header = T)
  
  ## Rename
  
  colnames(clone)<- newCols_v
  
  ## Add column with sample number
  snum <- (gsub(".*_S|*_alignment.*", '', file))
  clone$sample <- snum
  
  ## Return
  return(clone)
} # read_clone

### Read in clone files from multiple batches
readDir_clone<-function(dataDir, subDir, meta, norm_v) {
  #' Read in directory of files
  #' @description
  #' Given a base directory with multiple folders of clone files, read in all the files of a particular folder
  #' calls read_clone
  #' @param dataDir base directory containing multiple sub-directories of clone files
  #' @param subDir sub directory containing many clone files, can be all from one batch or from many batches
  #' @param meta metadata data.table used to identify batches if subDir contains clones from many batches
  #' @param norm_v character indicating if the raw clone count ("raw"), original norm counts ("orig"), or negative binomial norm count ("nb") should be used.
  #' @return List of lists of data.tables. One list for each batch in subDir, which contains 1 data.table for each file in that batch
  #' @export
  
  ## Get full directory path
  directory<-paste(dataDir,subDir, sep="")
  
  ## Get batches from metadata
  batchCol_v <- grep("atch", colnames(meta), value = T)
  batches_v <- as.character(unique(meta[,get(batchCol_v)]))
  
  ## Make empty list
  totalData_lslsdt <- list()
  
  ## Read in clones from each batch
  for (i in 1:length(batches_v)){
    ## Get batch
    currBatch_v <- batches_v[i]
    ## Get files
    currFiles_v <- list.files(directory, pattern = paste0(".*", currBatch_v, ".*txt$"))
    ## Sort by sample number
    currFiles_v <- currFiles_v[order(as.numeric(gsub("^.*_S|_align.*", "", currFiles_v)))]
    ## Get names
    currNames_v <- gsub(".*_S|_align.*", "", currFiles_v)
    ## Read in data
    #currData_lsdt <- lapply(currFiles_v, read_clone, directory = directory)
    currData_lsdt <- lapply(currFiles_v, function(x) read_clone(file = x, directory = directory, norm_v = norm_v))
    ## Add names
    names(currData_lsdt) <- currNames_v
    ## Add to list
    totalData_lslsdt[[currBatch_v]] <- currData_lsdt
  } # for i
  
  ## Return
  return(totalData_lslsdt)
} # readDir_clone

### Remove pseudo genes and fix V121/V122
cln_clones<-function(clonedf) {
  #' Clean Clone file
  #' @description
  #' Remove pseudo genes from clone data.table and also change V121 and V122 to be V1212
  #' @param clonedf data.table of a clonotype file
  #' @return data.table with pseudo genes removed 
  #' @export
  
  tbl<-clonedf
  pseudo=c("V22","V8","V10","V11","V123","V18","V21","V27","V28","V25")
  if (NROW(clonedf) != 0)
  {
    tbl<-tbl[!tbl$V %in% pseudo,]
    tbl[tbl=="V121"|tbl=="V122"]<-"V1212"  
  } 
  return(tbl)
}

### Merge list of data.tables
mergeDTs <- function(data_lsdt, mergeCol_v, keepCol_v) {
  #' Merge many data.tables together
  #' Take many data.tables and merge on and ID column, extracting a single column from each data.table as the column of interest
  #' @param data_lsdt list of data.tables to merge
  #' @param mergeCol_v which column from all of the data.tables to use to merge
  #' @param keepCol_v which column from all of the data.tables to use as the column of interest
  #' @value data.table with ncol == length(data_lsdt) + 1. Column names are names of list, or defaults to V1, V2,...
  #' @export
  
  ## Create initial table by extracting the 2 columns of interest from the rest
  merge_dt <- data_lsdt[[1]][,mget(c(mergeCol_v, keepCol_v))]
  
  ## Create initial column names (first check if list has names and add if not)
  if (is.null(names(data_lsdt))) {
    names_v	<- paste("V", 1:length(data_lsdt))
    names(data_lsdt) <- names_v
  } # fi
  
  colNames_v <- c(mergeCol_v, names(data_lsdt)[1])
  
  for	(i in 2:length(data_lsdt)) {
    ## Perform merge
    merge_dt <- merge(merge_dt,
                      data_lsdt[[i]][,mget(c(mergeCol_v, keepCol_v))],
                      by = mergeCol_v,
                      all = T)
    ## Update column names
    colNames_v <- c(colNames_v, names(data_lsdt)[i])
    ## Rename columns
    colnames(merge_dt) <- colNames_v
  } # for
  return(merge_dt)
} # mergeDTs

### Normalize 
norm4lib<-function(clonedf, median_v) {
  #' Perform median normalization
  #' @description
  #' Use batchwise median found separately to adjust all counts for a given sample
  #' @param clonedf data.table of all clonotypes for a sample
  #' @param median_v median value to adjust clonotype counts by
  #' @return data.table with normalized counts
  #' @export
  
  ## Get sum of all clones from a particular batch
  countSum_v<-sum(clonedf$CloneCount)
  
  ## Get scaling factor by dividing median by sum (will return how much greater or lesser this batch is than the median)
  scaleFactor_v<-median_v/countSum_v
  
  ## If no clones, return 0
  ## Otherwise, multiply raw count by factor
  if(countSum_v==0) {
    clonedf$normC=0
  } else {
    clonedf$normC<-round(clonedf$CloneCount * scaleFactor_v)
  } # fi
  
  ## Return
  return(clonedf)
}

### Group by V and J, summing count
groupVJ <- function(data_lsdt, column_v) {
  #' Group clone file by V and J
  #' @description
  #' Collapse all clonotype records for a file into 260 VJ combination records. Sum clone counts by shared VJ combination
  #' @param data_lsdt list of data.tables, each data.table contains all clonotype records for 1 sample
  #' @param column_v name of the count column to be aggregated
  #' @return list of data.tables grouped by VJ combinations
  #' @export
  
  lapply(data_lsdt, function(y) {
    ## Add vj column
    y$VJ <- paste(y$V, y$J, sep = "_")
    ## Remove other columns
    y[,(c("ID", "sample", "Seq", "aaSeq", "V", "J")) := NULL]
    ## Collapse by VJ column
    collapseY <- y[, sum(get(column_v)), by = .(VJ)]
    ## Return
    return(collapseY)
  })
} # groupVJ

### Get deviation between replicates in a batch
getDeviation <- function(compare_lsdt){
  #' Get deviation between replicates
  #' @description
  #' Takes a list of length 2, one element for each batch, and finds the deviation between replicates. NOTE - only works with 2 batches
  #' @param compare_lsdt list of length 2. Each element is a data.table with 260 rows (1 for each VJ combination)
  #' and 1 column for each sample from that batch. Column order in both data.tables corresponds to replicate identity (e.g. list[[1]][,1] is a replicate with list[[2]][,1])
  #' @return data.table with deviation ratios for each replicate pair
  #' @export
  
  ## Only works with 2 batches...
  batch1_dt <- compare_lsdt[[1]]
  batch2_dt <- compare_lsdt[[2]]
  ## Get replicates (plus one for VJ column)
  nReps_v <- unique(sapply(compare_lsdt, function(x) ncol(x)))
  ## Empty output data.table
  output_dt <- NULL
  
  for (i in 2:nReps_v){
    ## Merge Together
    currMerge_dt <- merge(batch1_dt[,c(1,i), with = F], batch2_dt[,c(1,i),with=F], by = "VJ", all = F)
    ## Colnames
    currCols_v <- colnames(currMerge_dt)
    ## Remove rows with 1 NA or value = 1
    naRows_v <- which( unlist(apply(currMerge_dt, 1, function(x) length(which(is.na(x))))) > 0)
    currMerge_dt <- currMerge_dt[-naRows_v,]
    ## Take logs
    currMerge_dt[, (currCols_v[-1]) := lapply(.SD, function(x) log2(x+1)), .SDcols = currCols_v[-1]]
    ## Ratio
    currMerge_dt$ratio <- currMerge_dt[,get(currCols_v[2])] - currMerge_dt[,get(currCols_v[3])]
    ## Add sample column
    currMerge_dt$pair <- paste(currCols_v[2], currCols_v[3], sep = "_")
    ## Rbind
    output_dt <- rbind(output_dt, currMerge_dt[,mget(c("VJ", "ratio", "pair"))])
  } # for i
  return(output_dt)
} # get Deviation
