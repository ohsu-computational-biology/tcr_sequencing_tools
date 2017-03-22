# QC Summary Script
# Aggregate all of the QC functionality into one script.
# 

# Load libraries
# .libPaths("/mnt/lustre1/CompBio/lib/R/library")
suppressMessages(source("https://bioconductor.org/biocLite.R", echo=FALSE, verbose=FALSE))
suppressMessages(library(ShortRead))
suppressMessages(library(stringr))
suppressMessages(library(Biostrings))

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

### Normalization
raw.clonotype.count.files <- arguments[9]
norm.clonotype.count.files <- arguments[10]
norm.output <- arguments[11]

######################################
### Count Spikes and Remove Spikes ###
######################################
# Creates QC report that was previously housed in the count.spikes script, and also creates
# the aggregate summary that was previously hosued in the count.spikes.QC script

# Separate input files
fastq.files <- strsplit(fastq.files, ',')
despiked.fastq.files <- strsplit(despiked.fastq.files, ',')
spike.count.files <- strsplit(spike.count.files, ',')

# Perform QC operations on each file

for (i in 1:length(fastq.files[[1]])){

    ### Read in data  
    fastq.reads <- readFastq(fastq.files[[1]][i])
    despiked.reads <- readFastq(despiked.fastq.files[[1]][i])
    spike.table <- read.table(spike.count.files[[1]][i], header = T, sep = ",")
    
    ### Calculate number of reads of each type
    num.fastqs <- length(fastq.reads)
    num.despiked <- length(despiked.reads)
    spiked.reads <- sum(spike.table$spike.count)

    ### Calculate various QC statistics
    
    # Total reads removed
    reads.removed <- num.fastqs - num.despiked
    
    # Percent reads retained
    pct.reads.retained <- round((num.despiked / num.fastqs * 100), digits = 1)
    
    # Percent spiked reads
    spike.percent <- (spiked.reads / num.fastqs) * 100
    
    # Aggregate into a single row data frame
    summary.data <- data.frame(list(fastq.files[[1]][i], num.fastqs, num.despiked,
                                    reads.removed, pct.reads.retained,
                                    spiked.reads, spike.percent))
    colnames(summary.data) <- c("ID", "Number Total Reads", "Number Reads After Despike",
                                "Reads Removed", "Percent original reads retained",
                                "Number Spiked Reads", "Percent Spiked Reads")
    
    # Extract spike counts and turn into single row data frame
    spike.out <- as.data.frame(spike.table$spike.count)
    rownames(spike.out) <- spike.table$SPIKE_ID
    spike.out <- t(spike.out)
    
    # Combine into one data frame
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
  
# Separate input files
align.reports <- strsplit(align.reports, ',')

# Create empty data frame
align.df <- data.frame()

# Summarize info from each report file
for (i in 1:length(align.reports[[1]])){
  
  # Extract information
  curr.record <- readLines(align.reports[[1]][i])
  curr.date <- paste(str_trim(str_split(curr.record[1], ":")[[1]][-1]), collapse='')
  curr.input <- basename(str_trim(str_split(curr.record[2], ':')[[1]][2]))
  curr.output <- basename(str_trim(str_split(curr.record[3], ':')[[1]][2]))
  curr.command <- str_trim(str_split(curr.record[4], ':')[[1]][2])
  curr.total <- str_trim(str_split(curr.record[5], ':')[[1]][2])
  curr.aligned <- str_trim(str_split(curr.record[6], ":")[[1]][2])
  subsequent <- NULL
  for (j in 7:14){
    curr.temp <- as.numeric(str_replace(str_trim(str_split(curr.record[j], ":")[[1]][2]),
                                        "%", ''))
    subsequent[j] <- curr.temp
  }
  subsequent <- subsequent[7:14]
  
  # Combine into a row to add to data frame
  row <- c(curr.date, curr.input, curr.output, curr.command, curr.total, curr.aligned,
           subsequent)

  # Populate data frame
  align.df <- rbind(align.df, row, stringsAsFactors = F)
  
}  #  for

# Add column names
  colnames(align.df) <- c("analysis.date", "inputs", "output", "command", "total.reads",
                          "aligned.reads", "aligned.pct", "filtered.diff.v.j.loci",
                          "failed.alignment.v.hits", "failed.alignment.j.hits", 
                          "failed.alignment.low.score", "overlapped.pct", 
                          "overlapped.and.aligned.pct", "overlapped.and.not.aligned.pct")

################
### Assembly ###
################
  
# Separate input files
assemble.reports <- strsplit(assemble.reports, ',')

# Create empty data frame
assemble.df <- data.frame()

# Summarize info from each report file
for (i in 1:length(assemble.reports[[1]])){
  
  # Extract information
  curr.record <- readLines(assemble.reports[[1]][i])
  curr.date <- paste(str_trim(str_split(curr.record[1], ":")[[1]][-1]), collapse='')
  curr.input <- basename(str_trim(str_split(curr.record[2], ':')[[1]][2]))
  curr.output <- basename(str_trim(str_split(curr.record[3], ':')[[1]][2]))
  curr.command <- str_trim(str_split(curr.record[4], ':')[[1]][2])
  curr.count <- str_trim(str_split(curr.record[5], ':')[[1]][2])
  curr.total <- str_trim(str_split(curr.record[6], ":")[[1]][2])
  subsequent <- NULL
  for (j in 7:14){
    curr.temp <- as.numeric(str_replace(str_trim(str_split(curr.record[j], ":")[[1]][2]),
                                        "%", ''))
    subsequent[j] <- curr.temp
  }
  subsequent <- subsequent[7:14]
  
  # Combine into a row to add to data frame
  row <- c(curr.date, curr.input, curr.output, curr.command, curr.count, curr.total,
           subsequent)
  
  # Populate data frame
  assemble.df <- rbind(assemble.df, row, stringsAsFactors = F)
  
}  #  for

# Add column names
colnames(assemble.df) <- c("analysis.date", "inputs", "output", "command", "clonotype.count",
                        "total.reads.used", "pct.total.reads.used", "pct.reads.used.as.core",
                        "pct.reads.used.mapped.low.quality", 
                        "pct.reads.used.PCR.error.corr.clustered", 
                        "clonotypes.elim.error.corr", "pct.reads.dropped.no.clonal.sequence", 
                        "pct.reads.dropped.low.quality", "pct.reads.dropped.failed.mapping")
  
  
#####################
### Normalization ###
#####################

# Separate input files
raw.clonotype.count.files <- strsplit(raw.clonotype.count.files, ',')
norm.clonotype.count.files <- strsplit(norm.clonotype.count.files, ',')

# Empty data frame
output.data <- NULL

for(i in 1:length(raw.clonotype.count.files[[1]]))  {
  #   Note that we use check.names=FALSE; this preserves the original column names,
  #       which is useful since some downstream tools (e.g. VDJTools' Convert() function)
  #       assume certain column names
  curr.raw <- read.table(raw.clonotype.count.files[[1]][i],
                         check.names=FALSE,  
                         stringsAsFactors=FALSE,
                         sep="\t",
                         header=TRUE);
  
  curr.normalized <- read.table(norm.clonotype.count.files[[1]][i],
                                check.names=FALSE,  
                                stringsAsFactors=FALSE,
                                sep="\t",
                                header=TRUE);
  
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
  #combined.table$raw.file.name <- raw.clone.counts[i];
  #combined.table$normalized.file.name <- processed.clone.counts[i];
  combined.table$normalization.factor <- round((combined.table$normalized.clone.count / combined.table$raw.clone.count), digits=1);
  
  #output.file.name <- paste(sample.id.raw.clone.counts[i], "_normalization_QC.txt", sep="");
  #output.file.name <- file.path(path.to.normalized.clone.counts, output.file.name);
  #cat("Writing output to: ", output.file.name, "\n", sep="");
  
  output.data <- rbind(output.data, summary(combined.table$normalization.factor))
  
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

### Normalization
write.table(output.data,
            file = norm.output,
            quote = FALSE,
            sep = '\t',
            row.names = FALSE)