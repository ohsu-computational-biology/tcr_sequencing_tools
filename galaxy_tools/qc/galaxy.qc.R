# QC Summary Script
# Aggregate all of the QC functionality into one script.
# 

# Load libraries
# .libPaths("/mnt/lustre1/CompBio/lib/R/library")
suppressMessages(source("https://bioconductor.org/biocLite.R", echo=FALSE, verbose=FALSE))
suppressMessages(library(ShortRead))
#suppressMessages(library(stringr))
suppressMessages(library(Biostrings))
suppressMessages(library(data.table))

#   disable scientific notation
options(scipen=999);

#################
### Arguments ###
#################

arguments <- commandArgs(trailingOnly=TRUE)

### Count Spikes and Remove Spikes
fastq.files <- arguments[1] # PEARed fastqs
spike.count.files <- arguments[2] # 25-bp spike counts
despiked.fastq.files <- arguments[3] # PEARed fastqs with spikes removed
spike.count.output <- arguments[4] # specifies output for QC table made by script

### MiXCR Alignment
align.reports <- arguments[5]
align.output <- arguments[6]

### MiXCR Assembly
assemble.reports <- arguments[7]
assemble.output <- arguments[8]

### Decontamination
decontam.qc.files <- arguments[9]
decontam.qc.out <- arguments[10]

### Normalization
decontam.clone.files <- arguments[11]
norm.clonotype.count.files <- arguments[12]
norm.output <- arguments[13]

######################################
### Count Spikes and Remove Spikes ###
######################################

### Creates QC report that was previously housed in the count.spikes script, and also creates
### the aggregate summary that was previously hosued in the count.spikes.QC script

### Separate input files
fastq.files <- unlist(strsplit(fastq.files, ','))
despiked.fastq.files <- unlist(strsplit(despiked.fastq.files, ','))
spike.count.files <- unlist(strsplit(spike.count.files, ','))

if (length(fastq.files) != length(despiked.fastq.files)){
    stop("Incorrect number of full and/or despiked files")
} else if (length(fastq.files) != length(spike.count.files)){
    stop("Incorrect number of full fastq and/or spike files")
} # fi


### Perform QC operations on each file

for (i in 1:length(fastq.files)){

    ## Read in data
    fastq.reads <- readFastq(fastq.files[i])
    despiked.reads <- readFastq(despiked.fastq.files[i])
    spike.table <- fread(spike.count.files[i])
    print(spike.table)
    
    ## Calculate number of reads of each type
    num.fastqs <- length(fastq.reads)
    num.despiked <- length(despiked.reads)
    spiked.reads <- sum(spike.table$spike.count)

    ## Calculate various QC statistics
    
    ## Total reads removed
    reads.removed <- num.fastqs - num.despiked
    
    ## Percent reads retained
    pct.reads.retained <- round((num.despiked / num.fastqs * 100), digits = 1)
    
    ## Percent spiked reads
    spike.percent <- (spiked.reads / num.fastqs) * 100

    cat("num.fastq\n", num.fastqs, '\n')
    cat("num.despiked\n", num.despiked, '\n')
    cat("spiked.reads\n", spiked.reads, "\n")
    cat("reads.removed", reads.removed, '\n')
    cat("pct.reads.retained", pct.reads.retained, '\n')
    cat("spike.percent", spike.percent, '\n')
    
    ## Aggregate into a single row data frame
    summary.data <- data.frame(list(fastq.files[i], num.fastqs, num.despiked,
                                    reads.removed, pct.reads.retained,
                                    spiked.reads, spike.percent))
    
    colnames(summary.data) <- c("ID", "Number Total Reads", "Number Reads After Despike",
                                "Reads Removed", "Percent original reads retained",
                                "Number Spiked Reads", "Percent Spiked Reads")
    
    ## Extract spike counts and turn into single row data frame
    spike.out <- as.data.frame(spike.table$spike.count)
    rownames(spike.out) <- spike.table$SPIKE_ID
    spike.out <- t(spike.out)
    
    ## Combine into one data frame
    new.row <- cbind(summary.data, spike.out)
    if ( i == 1){
        out.df <- new.row
    } else {
        out.df <- rbind(out.df, new.row)
    }   #   else
}   #   for
    
  
###################
### MiXCR Align ###
###################
  
### Separate input files
align.reports <- unlist(strsplit(align.reports, ','))

### Create empty data frame
align.df <- data.frame()

### Summarize info from each report file
for (i in 1:length(align.reports)){
  
    ## Extract information
    curr.record <- readLines(align.reports[i])
    curr.date <- paste(trimws(strsplit(curr.record[1], ":")[[1]][-1]), collapse='')
    curr.input <- basename(trimws(strsplit(curr.record[2], ':')[[1]][2]))
    curr.output <- basename(trimws(strsplit(curr.record[3], ':')[[1]][2]))
    curr.version <- trimws(strsplit(curr.record[4], ':|;')[[1]][2])
    curr.time <- trimws(strsplit(curr.record[5], ":")[[1]][2])
    curr.command <- gsub(" --report.*", '', trimws(strsplit(curr.record[6], ":")[[1]][2]))
    curr.total <- trimws(strsplit(curr.record[7], ':')[[1]][2])

    counter <- 0
    subsequent <- NULL
    for (j in 8:14){
        curr.temp <- unlist(strsplit(trimws(strsplit(curr.record[j], ":")[[1]][2]), ' '))
        curr.num <- as.numeric(curr.temp[1])
        curr.pct <- as.numeric(gsub("\\(|%|\\)", "", curr.temp[2]))
        subsequent[j+counter] <- curr.num
        subsequent[j+counter+1] <- curr.pct
        counter <- counter + 1
    } # for j

    ## Remove leading NAs
    subsequent <- subsequent[8:21]
  
  ## Combine into a row to add to data frame
  row <- c(curr.date, curr.input, curr.output, curr.version, curr.time, curr.command,
           curr.total, subsequent)

  # Populate data frame
  align.df <- rbind(align.df, row, stringsAsFactors = F)
  
}  #  for

### Add column names
colnames(align.df) <- c("analysis.date", "inputs", "output", "version", "time",
                        "command", "total.reads", "aligned.reads", "aligned.pct", "failed.alignment.no.hits",
                        "pct.no.hits", "failed.alignment.no.j.hits", "pct.no.j.hits", "failed.alignment.low.score",
                        "pct.low.score", "num.overlapped", "pct.overlapped", "num.overlapped.and.aligned",
                        "pct.overlapped.and.aligned", "num.overlapped.and.not.aligned", "pct.overlapped.and.not.aligned")

################
### Assembly ###
################
  
### Separate input files
assemble.reports <- unlist(strsplit(assemble.reports, ','))

### Create empty data frame
assemble.df <- data.frame()

### Summarize info from each report file
for (i in 1:length(assemble.reports)){
  
    ## Extract information
    curr.record <- readLines(assemble.reports[i])
    curr.date <- paste(trimws(strsplit(curr.record[1], ":")[[1]][-1]), collapse='')
    curr.input <- basename(trimws(strsplit(curr.record[2], ':')[[1]][2]))
    curr.output <- basename(trimws(strsplit(curr.record[3], ':')[[1]][2]))
    curr.version <- trimws(strsplit(curr.record[4], ":|;")[[1]][2])
    curr.time <- trimws(strsplit(curr.record[5], ":")[[1]][2])
    curr.command <- trimws(strsplit(curr.record[6], ':')[[1]][2])
    curr.count <- trimws(strsplit(curr.record[7], ':')[[1]][2])
    curr.avg.per.clone <- trimws(strsplit(curr.record[8], ":")[[1]][2])

    subsequent <- NULL
    for (j in 9:18){
        curr.temp <- unlist(strsplit(trimws(strsplit(curr.record[j], ":")[[1]][2]), " "))
        curr.num <- as.numeric(curr.temp[1])
        curr.pct <- as.numeric(gsub("\\(|%|\\)", "", curr.temp[2]))
        subsequent[j+counter] <- curr.num
        subsequent[j+counter+1] <- curr.pct
        counter <- counter + 1
    } # for j

    for (j in 19:21){
        curr.temp <- as.numeric(trimws(strsplit(curr.record[j], ":")[[1]][2]))
        subsequent[j+counter] <- curr.temp
    } # for j

    ## Remove NAs
    subsequent <- subsequent[9:31]
  
  # Combine into a row to add to data frame
    row <- c(curr.date, curr.input, curr.output, curr.version, curr.time, curr.command,
             curr.count, curr.avg.per.clone, subsequent)
  
  # Populate data frame
  assemble.df <- rbind(assemble.df, row, stringsAsFactors = F)
  
}  #  for

# Add column names
colnames(assemble.df) <- c("analysis.date", "inputs", "output", "version", "time", "command", "clonotype.count",
                           "avg.reads.per.clonotype", "num.reads.used", "pct.used.of.total", "num.reads.used.b4.clust",
                           "pct.of.total", "num.reads.used.as.core", "pct.of.used", "num.reads.mapped.low.quality", "pct.mapped.of.used",
                           "num.PCR.error.clust", "pct.PCR.clust.of.used", "num.VJC.clust", "pct.VJC.clust.of.used",
                           "num.dropped.no.clonal.seq", "pct.dropped.no.clonal", "num.reads.dropped.low.qual",
                           "pct.dropped.low.quality", "num.reads.dropped.fail.map", "pct.dropped.fail.map",
                           "num.reads.dropped.low.qual.clone", "pct.dropped.low.qual.clone", "clonotypes.elim.error.corr",
                           "clonotypes.dropped.low.qual", "clonotypes.pre.clust.similar.VJC")

#######################
### Decontamination ###
#######################

### Get files
decontam_qc_files_v <- unlist(strsplit(decontam.qc.files, ','))

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

# Separate input files
raw.clonotype.count.files <- unlist(strsplit(decontam.clone.files, ','))
norm.clonotype.count.files <- unlist(strsplit(norm.clonotype.count.files, ','))

# Empty data frame
output.data <- NULL

for(i in 1:length(raw.clonotype.count.files))  {
  #   Note that we use check.names=FALSE; this preserves the original column names,
  #       which is useful since some downstream tools (e.g. VDJTools' Convert() function)
  #       assume certain column names
  curr.raw <- fread(raw.clonotype.count.files[i])
  
  curr.normalized <- fread(norm.clonotype.count.files[i])
  
  # Attempt at parallelization test
  #mismatches <- which(curr.raw$`Clonal sequence(s)` != curr.normalized$`Clonal sequence(s)`)
  #if (length(mismatches) == 0) {
  #  stop("Clonal sequences do not match. Check parallelism of files.")
  #}  #  if
  
  #   Basic QC
  curr.raw.CDR3 <- curr.raw$`AA. Seq. CDR3`;
  curr.normalized.CDR3 <- curr.normalized$"AA. Seq. CDR3";
  if(!identical(curr.raw.CDR3, curr.normalized.CDR3)) {
    stop("Mistmatch between amino acid CDR3 region, raw and normalized");
  }   #   fi
  
  combined.table <- data.frame(raw.clone.count=curr.raw$"Clone count");
  combined.table$normalized.clone.count <- curr.normalized$"Normalized clone count";
  combined.table$raw.clone.percent <- curr.raw$"Clone fraction";
  combined.table$normalized.clone.percent <- curr.normalized$"Normalized clone fraction";
  combined.table$normalization.factor <- round((combined.table$normalized.clone.count / combined.table$raw.clone.count), digits=1);
  
  #output.file.name <- paste(sample.id.raw.clone.counts[i], "_normalization_QC.txt", sep="");
  #output.file.name <- file.path(path.to.normalized.clone.counts, output.file.name);
  #cat("Writing output to: ", output.file.name, "\n", sep="");
  
    output.data <- rbind(output.data, c(apply(combined.table[,1:4], 2, mean) , as.vector(summary(combined.table$normalization.factor))))
    colnames(output.data) <- c("avg.raw.count", "avg.norm.count", "avg.raw.pct", "avg.norm.pct", "min.scale", "lower.q.scale",
                               "median.scale", "mean.scale", "upper.q.scale", "max.scale")
  
  #write.table(combined.table,
  #            file=output.file.name,
  #            quote=FALSE,
  #            sep=",",
  #            row.names=FALSE);
  
  #   reset value
  rm(combined.table);
}   #   for i

###############
### Outputs ###
###############

### Spike Count & Remove Spikes
write.table(out.df, 
            file=spike.count.output,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

### MiXCR Align
write.table(align.df,
            file = align.output,
            quote = FALSE,
            sep = '\t',
            row.names = FALSE)

### MiXCR Assemble
write.table(assemble.df,
            file = assemble.output,
            quote = FALSE,
            sep = '\t',
            row.names = FALSE)

### Decontamination
write.table(decontam_output_matrix,
            file = decontam.qc.out,
            sep = '\t',
            quote = FALSE,
            row.names = F)

### Normalization
write.table(output.data,
            file = norm.output,
            quote = FALSE,
            sep = '\t',
            row.names = FALSE)
