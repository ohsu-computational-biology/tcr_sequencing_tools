###
### QC Summary Script
###

### Aggregate all of the QC functionality into one script.
  ### PEAR
  ### Count Spikes - 25-bp
  ### Spike Removal
  ### MiXCR
  ### Decontamination
  ### Normalization

### disable scientific notation
options(scipen=999);
suppressMessages(library(data.table))
suppressWarnings(suppressMessages(library(ShortRead)))
suppressWarnings(suppressMessages(library(Biostrings)))

#################
### Arguments ###
#################

arguments <- commandArgs(trailingOnly=TRUE)

### Count Spikes and Remove Spikes
fastqFiles_v <- arguments[1]
fastqNames_v <- arguments[2]

spikeCountFiles_v <- arguments[3]
spikeCountNames_v <- arguments[4]

despikedFiles_v <- arguments[5]
despikedNames_v <- arguments[6]

spikeCountOut_v <- arguments[7]

### MiXCR Align and Assemble
alignReports_v <- arguments[8]
alignNames_v <- arguments[9]
alignOut_v <- arguments[10]

assembleReports_v <- arguments[11]
assembleNames_v <- arguments[12]
assembleOut_v <- arguments[13]

### Decontamination
decontamQCFiles_v <- arguments[14]
decontamQCNames_v <- arguments[15]
decontamOut_v <- arguments[16]

### Normalization
decontamCloneFiles_v <- arguments[17]
decontamCloneNames_v <- arguments[18]

normCountFiles_v <- arguments[19]
normCountNames_v <- arguments[20]

normOut_v <- arguments[21]

debug_v <- arguments[22]

###################
### For Testing ###
###################

# fastqFiles_v <- "~/galaxy/test-data/count.spikes/pear_S10_data.fastq"
# fastqNames_v <- "~/galaxy/test-data/qc/temp/fastqNames.txt"
# 
# spikeCountFiles_v <- "~/galaxy/test-data/count.spikes/25bp_S10_counts.txt"
# spikeCountNames_v <- "~/galaxy/test-data/qc/temp/spikeCountNames.txt"
# 
# despikedFiles_v <- "~/galaxy/test-data/remove.spikes/reads_S10_data.fastq"
# despikedNames_v <- "~/galaxy/test-data/qc/temp/despikedNames.txt"
# 
# spikeCountOut_v <- "testSpikeCountOut.txt"
# 
# alignReports_v <- "~/galaxy/test-data/mixcr/mixcr_S10_alignReport.txt"
# alignNames_v <- "~/galaxy/test-data/qc/temp/alignNames.txt"
# alignOut_v <- "testAlignOut.txt"
# 
# assembleReports_v <- "~/galaxy/test-data/mixcr/mixcr_S10_assembleReport.txt"
# assembleNames_v <- "~/galaxy/test-data/qc/temp/assembleNames.txt"
# assembleOut_v <- "testAssembleOut.txt"
# 
# decontamQCFiles_v <- "~/galaxy/test-data/decontam/decontam_S10_qc.txt"
# decontamQCNames_v <- "~/galaxy/test-data/qc/temp/decontamQCNames.txt"
# decontamOut_v <- "testDecontamOut.txt"
# 
# decontamCloneFiles_v <- "~/galaxy/test-data/decontam/decontam_S10_clones.txt"
# decontamCloneNames_v <- "~/galaxy/test-data/qc/temp/decontamCloneNames.txt"
# 
# normCountFiles_v <- "~/galaxy/test-data/normalization/normalize_S10_clones.txt"
# normCountNames_v <- "~/galaxy/test-data/qc/temp/normNames.txt"
# 
# normOut_v <- "testNormOut.txt"

#################
### Functions ###
#################
dsName_v <- function(dataSetName_v){
  # Return last directory from file path and file name of galaxy dataset.dat path
  # dataSetName_v - single element of a galaxy dataset path (e.g. /database/files/000/dataset123.dat)
  splitName_v <- strsplit(dataSetName_v, split = "/")[[1]]
  x_v <- length(splitName_v)
  newName_v <- paste(splitName_v[c((x_v-1),x_v)], collapse = "/")
  return(newName_v)
}

extractRecord <- function(string_v, record_v, names_v){
  # Extract the number and percent values from mixcr records
  # string_v - character vector to distinguish which element of record to look at
  # record_v - assemblage of character records read in from mixcr record
  # names_v - character vectors to use as column names
  grep_v <- grep(string_v, record_v, value = T)
  if (length(grep_v) == 0){
    final_v <- c(NA, NA)
    names(final_v) <- names_v
  } else {
    partial_v <- unlist(strsplit(trimws(strsplit(grep_v, ":")[[1]][2]), " "))
    final_v <- as.numeric(gsub("\\(|%|\\)", "", partial_v))
    names(final_v) <- names_v
  }
  return(final_v)
}

######################################
### Count Spikes and Remove Spikes ###
######################################

### Creates QC report that was previously housed in the count.spikes script, and also creates
### the aggregate summary that was previously hosued in the count.spikes.QC script

### Separate input files
fastqFiles_v <- unlist(strsplit(fastqFiles_v, ','))
despikedFiles_v <- unlist(strsplit(despikedFiles_v, ','))
spikeCountFiles_v <- unlist(strsplit(spikeCountFiles_v, ','))

### Read in names
fastqNames_dt <- fread(fastqNames_v, header = F);
spikeCountNames_dt <- fread(spikeCountNames_v, header = F);
despikedNames_dt <- fread(despikedNames_v, header = F); 

### Check Names
if (debug_v){
  print("Fastq Names: ")
  print(fastqNames_dt)
  print("Spike Count Names: ")
  print(spikeCountNames_dt)
  print("Despiked Fastq Names: ")
  print(despikedNames_dt)
}

### Compare
fastqToSpikeCompare_v <- which(fastqNames_dt$V1 != spikeCountNames_dt$V1)
spikeToDespikeCompare_v <- which(spikeCountNames_dt$V1 != despikedNames_dt$V1)

if (length(fastqToSpikeCompare_v) > 0 | length(spikeToDespikeCompare_v) > 0){
  stop("PEAR fastq, spike count, and despiked fastq files do not match")
}

### Perform QC operations on each file
if (debug_v) print("Begin Spike Count/Remove")

for (i in 1:length(fastqFiles_v)){
  if (debug_v) {print("Currently on: "); print(i)}

    ## Read in data
    fastqReads_srq <- readFastq(fastqFiles_v[i])
    despikedReads_srq <- readFastq(despikedFiles_v[i])
    spikeTable_dt <- fread(spikeCountFiles_v[i])
    if (debug_v) {print("Spike Table: "); print(head(spikeTable_dt))}

    ## Get File information
    currFastqDatasetName_v <- dsName_v(fastqFiles_v[i])
    currDespikedDatasetName_v <- dsName_v(despikedFiles_v[i])
    currSpikeCountDatasetName_v <- dsName_v(spikeCountFiles_v[i])
    
    ## Get sample names
    currFastqSampleName_v <- fastqNames_dt$V1[i]
        
    ## Calculate number of reads of each type
    numFastqs_v <- length(fastqReads_srq)
    numDespiked_v <- length(despikedReads_srq)
    spikedReads_v <- sum(spikeTable_dt$spike.count)

    ##
    ## Calculate various QC statistics
    ##

    ## Total reads removed
    numRemoved_v <- numFastqs_v - numDespiked_v

    ## Percent reads retained
    pctReadsKept_v <- round((numDespiked_v / numFastqs_v * 100), digits = 1)

    ## Percent spiked reads
    pctSpiked_v <- (spikedReads_v / numFastqs_v) * 100
    
    if (debug_v){
      cat("num.fastq\n", numFastqs_v, '\n')
      cat("num.despiked\n", numDespiked_v, '\n')
      cat("spiked.reads\n", spikedReads_v, "\n")
      cat("reads.removed", numRemoved_v, '\n')
      cat("pct.reads.retained", pctReadsKept_v, '\n')
      cat("spike.percent", pctSpiked_v, '\n')
      
      print("Begin df")
    }
    ## Aggregate into a single row data frame
    summaryData_df <- data.frame("Sample" = currFastqSampleName_v, "Total_Reads" = numFastqs_v,
                                 "Despiked_Reads" = numDespiked_v, "Reads_Removed" = numRemoved_v,
                                 "Pct_Kept" = pctReadsKept_v, "Spiked_Reads" = spikedReads_v,
                                 "Pct_Spike" = pctSpiked_v, "Galaxy_Fastq_Name" = currFastqDatasetName_v,
                                 "Galaxy_Despiked_Name" = currDespikedDatasetName_v, 
                                 "Galaxy_SpikeCountName_v" = currSpikeCountDatasetName_v)
    if (debug_v) print("end summary df")
    
    ## Extract spike counts and turn into single row data frame
    spikeOut_df <- as.data.frame(spikeTable_dt$spike.count)
    rownames(spikeOut_df) <- spikeTable_dt$SPIKE_ID
    spikeOut_df <- t(spikeOut_df)
    
    ## Combine into one data frame
    newRow_df <- cbind(summaryData_df, spikeOut_df)
    if ( i == 1){
      if (debug_v) print("i == 1")
        summaryOut_df <- newRow_df
        if (debug_v) {print(summaryOut_df); print("Done i == 1")}
    } else {
      if (debug_v) print("i != 1")
        summaryOut_df <- rbind(summaryOut_df, newRow_df)
        if (debug_v) print("Done i != 1")
    }   #   else
    
}   #   for

###################
### MiXCR Align ###
###################

if (debug_v) print("Begin MiXCR Align")

### Separate input files
alignReports_v <- unlist(strsplit(alignReports_v, ','))

### Create empty data frame
align.df <- data.frame()

### Summarize info from each report file
for (i in 1:length(alignReports_v)){

    ## Get Record
    curr.record <- readLines(alignReports_v[i])
    
    ## Date
    curr.date <- grep("Analysis Date", curr.record, value = T)
    curr.date <- trimws(gsub("[A-z ]+:|[0-9]+:.*PDT ", "", curr.date))
    names(curr.date) <- "analysis.date"
    
    ## I/0
    curr.input <- basename(trimws(strsplit(grep("Input file", curr.record, value = T), ':')[[1]][2]))
    curr.output <- basename(trimws(strsplit(grep("Output file", curr.record, value = T), ':')[[1]][2]))
    names(curr.input) <- "inputs"; names(curr.output) <- "output"
    
    ## Version
    curr.version <- trimws(strsplit(grep("Version", curr.record, value = T), ":")[[1]][2])
    names(curr.version) <- "version"
    
    ## Total Reads
    curr.total <- trimws(strsplit(grep("Total seq", curr.record, value = T), ':')[[1]][2])
    names(curr.total) <- "total.reads"

    ## Alignment
    curr.Success <- extractRecord("Success", curr.record, c("aligned.reads", "aligned.pct"))
    curr.NoHit <- extractRecord("no hits", curr.record, c("failed.align.no.hits", "pct.no.hits"))
    curr.NoJ <- extractRecord("J hits", curr.record, c("failed.align.no.j", "pct.no.j"))
    curr.Low <- extractRecord("low total", curr.record, c("failed.aling.low.score", "pct.low.score"))
    curr.Overlap <- extractRecord("Overlapped:", curr.record, c("num.overlapped", "pct.overlapped"))
    curr.Overlap.Align <- extractRecord("Overlapped and aligned:", curr.record,
                                        c("num.overlapped.and.aligned", "pct.overlapped.and.aligned"))
    curr.Overlap.Not.Align <- extractRecord("Overlapped and not", curr.record,
                                            c("num.overlapped.and.not.algned", "pct.overlapped.and.not.aligned"))

  ## Combine into a row to add to data frame
    row_v <- c(curr.date, curr.input, curr.output, curr.version, curr.total,
               curr.Success, curr.NoHit, curr.NoJ, curr.Overlap, curr.Overlap.Align, curr.Overlap.Not.Align)
    colNames_v <- names(row_v)
  
  # Populate data frame
  align.df <- rbind(align.df, row_v, stringsAsFactors = F)
  colnames(align.df) <- colNames_v

}  #  for

################
### Assembly ###
################

if (debug_v) print("Begin MiXCR Assemble")

### Separate input files
assembleReports_v <- unlist(strsplit(assembleReports_v, ','))

### Create empty data frame
assemble.df <- data.frame()

### Summarize info from each report file
for (i in 1:length(assembleReports_v)){

    ## Get Record
    curr.record <- readLines(assembleReports_v[i])
    
    ## Date
    curr.date <- grep("Analysis Date", curr.record, value = T)
    curr.date <- trimws(gsub("[A-z ]+:|[0-9]+:.*PDT ", "", curr.date))
    names(curr.date) <- "analysis.date"
    
    ## I/0
    curr.input <- basename(trimws(strsplit(grep("Input file", curr.record, value = T), ':')[[1]][2]))
    curr.output <- basename(trimws(strsplit(grep("Output file", curr.record, value = T), ':')[[1]][2]))
    names(curr.input) <- "inputs"; names(curr.output) <- "output"
    
    ## Version
    curr.version <- trimws(strsplit(grep("Version", curr.record, value = T), ":")[[1]][2])
    names(curr.version) <- "version"
    
    ## Count
    curr.count <- trimws(strsplit(grep("Final clonotype count", curr.record, value = T), ":")[[1]][2])
    curr.avg.per.clone <- trimws(strsplit(grep("Average number", curr.record, value = T), ":")[[1]][2])
    names(curr.count) <- "clonotype.count"; names(curr.avg.per.clone) <- "avg.reads.per.clonotype"
    
    curr.reads.used <- extractRecord("clonotypes, percent", curr.record, c("num.reads.used", "pct.used.of.total"))
    curr.reads.cluster <- extractRecord("clonotypes before", curr.record, c("num.reads.used.b4.clust", "pct.of.total"))
    curr.core <- extractRecord("used as a core", curr.record, c("num.reads.used.as.core", "pct.of.used"))
    curr.low <- extractRecord("quality reads", curr.record, c("num.reads.mapped.lowq", "pct.mapped.of.used"))
    curr.clust <- extractRecord("Reads clustered", curr.record, c("num.PCR.error.clust", "pct.PCR.clust.of.used"))
    curr.pre.clust <- extractRecord("pre-clustered", curr.record, c("num.VJC.clust", "pct.VJC.clust.of.used"))
    curr.dropped.lack <- extractRecord("lack of a clone", curr.record, c("num.drop.no.clonal.seq", "pct.dropped.no.clonal"))
    curr.dropped.low <- extractRecord("dropped due to low", curr.record, c("num.drop.lowq", "pct.dropped.lowq"))
    curr.dropped.fail <- extractRecord("failed mapping", curr.record, c("num.drop.fail.map", "pct.dropped.fail.map"))
    curr.dropped.low.clone <- extractRecord("low quality clones", curr.record, c("num.drop.lowq.clone", "pct.dropped.lowq.clone"))
    curr.pcr.correct <-  extractRecord("eliminated by", curr.record, "clonotypes.elim.PCR.error")
    curr.clone.dropped.lowq <- extractRecord("Clonotypes dropped", curr.record, "clonotypes.drop.lowq")
    curr.clone.preclust <- extractRecord("Clonotypes pre-clustered", curr.record, "clonotypes.pre.clust.similar.VJC")
    

  # Combine into a row to add to data frame
    row_v <- c(curr.date, curr.input, curr.output, curr.version, curr.count, curr.avg.per.clone, 
               curr.reads.used, curr.reads.cluster, curr.core, curr.low, curr.clust, curr.pre.clust,
               curr.dropped.lack, curr.dropped.low, curr.dropped.fail, curr.dropped.low.clone,
               curr.pcr.correct, curr.clone.dropped.lowq, curr.clone.preclust)
    colNames_v <- names(row_v)

  # Populate data frame
  assemble.df <- rbind(assemble.df, row_v, stringsAsFactors = F)
  colnames(assemble.df) <- colNames_v

}  #  for

#######################
### Decontamination ###
#######################

if (debug_v) print("Begin Decontam")

### Get files
decontam_qc_files_v <- unlist(strsplit(decontamQCFiles_v, ','))

decontam_output_matrix <- matrix(nrow = length(decontam_qc_files_v), ncol = ncol(fread(decontam_qc_files_v[1], nrows = 0)))

for (i in 1:length(decontam_qc_files_v)){
    ## Read data
    curr_data_dt <- fread(decontam_qc_files_v[i])
    ## Add to matrix
    decontam_output_matrix[i,] <- unlist(curr_data_dt[1,], use.names = F)
} # for

colnames(decontam_output_matrix) <- colnames(curr_data_dt)

#####################
### Normalization ###
#####################

if (debug_v) print("Begin Normalization")

# Separate input files
raw.clonotype.count.files <- unlist(strsplit(decontamCloneFiles_v, ','))
normCountFiles_v <- unlist(strsplit(normCountFiles_v, ','))

# Empty data frame
output.data <- NULL

for(i in 1:length(raw.clonotype.count.files))  {
  ## Read in raw and norm
  curr.raw <- fread(raw.clonotype.count.files[i])
  curr.normalized <- fread(normCountFiles_v[i])

  ## Get column names
  cdr3Col_v <- grep("AA. Seq. CDR3|aaSeqCDR3", colnames(curr.raw), value = T)
  rawCountCol_v <- grep("Clone count|cloneCount", colnames(curr.raw), value = T)
  rawFreqCol_v <- grep("Clone fraction|cloneFraction", colnames(curr.raw), value = T)
  normCountCol_v <- grep("Normalized clone count", colnames(curr.normalized), value = T)
  normFreqCol_v <- grep("Normalized clone fraction", colnames(curr.normalized), value = T)
  nbCountCol_v <- grep("nb.clone.count", colnames(curr.normalized), value = T)
  nbFreqCol_v <- grep("nb.clone.fraction", colnames(curr.normalized), value = T)
  
  ## Check if files are same
  curr.raw.CDR3 <- curr.raw[[cdr3Col_v]]
  curr.normalized.CDR3 <- curr.normalized[[cdr3Col_v]]
  if(!identical(curr.raw.CDR3, curr.normalized.CDR3)) {
    stop("Mistmatch between amino acid CDR3 region, raw and normalized");
  }   #   fi

  ## QC Table
  combined.table <- data.frame(raw.clone.count=curr.raw[[rawCountCol_v]]);
  ncols_v <- 2
  if (length(normCountCol_v) > 0){
    combined.table$median.norm.count <- curr.normalized[[normCountCol_v]]
    ncols_v <- ncols_v + 2
  } # fi
  if (length(nbCountCol_v) > 0){
    combined.table$nb.clone.count <- curr.normalized[[nbCountCol_v]]
    ncols_v <- ncols_v + 2
  }
  combined.table$raw.clone.percent <- curr.raw[[rawFreqCol_v]]
  if (length(normFreqCol_v) > 0){
    combined.table$median.norm.fraction <- curr.normalized[[normFreqCol_v]]
  } # fi
  if (length(nbFreqCol_v) > 0){
    combined.table$nb.clone.fraction <- curr.normalized[[nbFreqCol_v]]
  }
  if (length(normFreqCol_v) > 0){
    combined.table$median.norm.factor <- round((combined.table$median.norm.count / combined.table$raw.clone.count), digits = 1)
  }
  if (length(nbFreqCol_v) > 0){
    combined.table$nb.norm.factor <- round((combined.table$nb.clone.count / combined.table$raw.clone.count), digits = 1)
  }
  
  ## Summarize
  row_v <- apply(combined.table[,1:ncols_v], 2, mean)
  if (length(normFreqCol_v) > 0) {
    summary_v <- as.vector(summary(combined.table$median.norm.factor))[c(1,3,6)]
    names(summary_v) <- c("min.scale.median.norm", "median.scale.median.norm", "max.scale.median.norm")
    row_v <- c(row_v, summary_v)
  }
  if (length(nbFreqCol_v) > 0) {
    summary_v <- as.vector(summary(combined.table$nb.norm.factor))[c(1,3,6)]
    names(summary_v) <- c("min.scale.nb.norm", "median.scale.nb.norm", "max.scale.nb.norm")
    row_v <- c(row_v, summary_v)
  }
  
  ## Combine
  colNames_v <- names(row_v)
  output.data <- rbind(output.data, row_v)
  colnames(output.data) <- colNames_v

  rm(combined.table);
}   #   for i

###############
### Outputs ###
###############

if (debug_v) print("Begin Output")

### Spike Count & Remove Spikes
write.table(summaryOut_df,
            file=spikeCountOut_v,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

### MiXCR Align
write.table(align.df,
            file = alignOut_v,
            quote = FALSE,
            sep = '\t',
            row.names = FALSE)

### MiXCR Assemble
write.table(assemble.df,
            file = assembleOut_v,
            quote = FALSE,
            sep = '\t',
            row.names = FALSE)

### Decontamination
write.table(decontam_output_matrix,
            file = decontamOut_v,
            sep = '\t',
            quote = FALSE,
            row.names = F)

### Normalization
write.table(output.data,
            file = normOut_v,
            quote = FALSE,
            sep = '\t',
            row.names = FALSE)
