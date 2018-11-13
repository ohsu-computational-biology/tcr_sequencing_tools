#!/usr/bin/Rscript
### RScript that aggregates reports from MiXCR's assembly tool
### Basically turn their report files into a single text file that can be easily manipulated.

####################
### DEPENDENCIES ###
####################

library(optparse)
source(file.path(Sys.getenv("tool"), "misc/helperFxn.R"))
options(warn=1)

####################
### COMMAND ARGS ###
####################

optlist <- list(
        make_option(
                c("-i", "--inputDir"),
                type = "character",
                help = "Path to directory containing all and only MiXCR alignment reports for a batch."
        ),
        make_option(
                c("-o", "--outDir"),
                type = "character",
                help = "Path to directory where aggregated results will be written. File is mixcr.alignment.QC.summary.txt"
        )
)

##################
### PARSE ARGS ###
##################

p <- OptionParser(usage = "%prog -i inputDir -o outDir",
                option_list = optlist)
args <- parse_args(p)
opt <- args$options

inputDir_v <- args$inputDir
outDir_v <- args$outDir

###############
### ACTIONS ###
###############

### Get files
inputFiles_v <- list.files(inputDir_v)

### Sort files
inputFiles_v <- inputFiles_v[order(as.numeric(gsub("^S|_assemble.*", "", inputFiles_v)))]

### Create empty data.frame
output_df <- data.frame()

### Summarize info from each report file
for (i in 1:length(inputFiles_v)){

    ## Get Record
    currFile_v <- inputFiles_v[i]
    curr.record <- readLines(file.path(inputDir_v, currFile_v))
    
    ## Update user
    if (length(curr.record) != 22){
	print(sprintf("Sample %s has %d lines. The norm is 22.", currFile_v, length(curr.record)))
    }

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
  output_df <- rbind(output_df, row_v, stringsAsFactors = F)
  colnames(output_df) <- colNames_v

}  #  for

### Write
write.table(output_df, 
            file=file.path(outDir_v, "mixcr.assemble.QC.summary.txt"),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

