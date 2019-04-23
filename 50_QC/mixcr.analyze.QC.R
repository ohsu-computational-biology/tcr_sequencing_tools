#!/usr/bin/Rscript
### RScript that aggregates reports from MiXCR's alignment tool
### Basically turn their report files into a single text file that can be easily manipulated

####################
### DEPENDENCIES ###
####################

library(optparse)
source(file.path(Sys.getenv("tool"), "50_QC/helperFxn.R"))
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
inputFiles_v <- inputFiles_v[order(as.numeric(gsub("^S|_report.*", "", inputFiles_v)))]

### Create empty data.frame
align_df <- assemble_df <- data.frame()

### Extract info from each report
for (i in 1:length(inputFiles_v)){

    ## Get Record
    currFile_v <- inputFiles_v[i]
    curr.record <- readLines(file.path(inputDir_v, currFile_v))

    ## Update user
    if (!(length(curr.record) %in% c(39,40,41))){
	print(sprintf("Sample %s has %d lines. The norm is 39-41.", currFile_v, length(curr.record)))
    }

    ## Date
    curr.date <- grep("Analysis date", curr.record, value = T)[1]
    curr.date <- trimws(gsub("[A-z ]+:|[0-9]+:.*PDT ", "", curr.date))
    names(curr.date) <- "analysis.date"
    
    ## I/0
    curr.input <- basename(trimws(strsplit(grep("Input file", curr.record, value = T)[1], ':')[[1]][2]))
    curr.output <- basename(trimws(strsplit(grep("Output file", curr.record, value = T)[1], ':')[[1]][2]))
    curr.output2 <- basename(trimws(strsplit(grep("Output file", curr.record, value = T)[2], ':')[[1]][2]))
    names(curr.input) <- "inputs"; names(curr.output) <- "vdjca"; names(curr.output2) <- "clna"
    
    ## Version
    curr.version <- trimws(strsplit(grep("Version", curr.record, value = T)[1], ":")[[1]][2])
    names(curr.version) <- "version"
    
    ## Total Reads
    curr.total <- trimws(strsplit(grep("Total seq", curr.record, value = T), ':')[[1]][2])
    names(curr.total) <- "total.reads"

    ## Alignment
    curr.Success <- extractRecord("Success", curr.record, c("aligned.reads", "aligned.pct"))
    curr.NoHit <- extractRecord("no hits", curr.record, c("failed.align.no.hits", "pct.no.hits"))
    curr.NoJ <- extractRecord("J hits", curr.record, c("failed.align.no.j", "pct.no.j"))
    curr.Low <- extractRecord("low total", curr.record, c("failed.align.low.score", "pct.low.score"))
    curr.Overlap <- extractRecord("Overlapped:", curr.record, c("num.overlapped", "pct.overlapped"))
    curr.Overlap.Align <- extractRecord("Overlapped and aligned:", curr.record,
                                        c("num.overlapped.and.aligned", "pct.overlapped.and.aligned"))
    curr.Alignment.Aided <- extractRecord("Alignment-aided", curr.record, 
					c("num.align.aided.overlap", "pct.align.aided.overlap"))
    curr.Overlap.Not.Align <- extractRecord("Overlapped and not", curr.record,
                                            c("num.overlapped.and.not.algned", "pct.overlapped.and.not.aligned"))
    curr.TRB.chains <- extractRecord("TRB chains", curr.record, c("num.TRB.chains", "pct.TRB.chains"))[1]

    ## Assembly
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
    curr.clone.TRB.chains <- extractRecord("TRB chains", curr.record, c("num.clone.TRB.chains", "pct.clone.TRB.chains"))[2]

  ## Combine into a row to add to data frame
    alignRow_v <- c(curr.date, curr.input, curr.output, curr.output2, curr.version, curr.total,
               curr.Success, curr.NoHit, curr.NoJ, curr.Overlap, curr.Overlap.Align, curr.Overlap.Not.Align)
    alignNames_v <- names(alignRow_v)

    assembleRow_v <- c(curr.date, curr.input, curr.output,, curr.output2, curr.version, curr.count, curr.avg.per.clone,
               curr.reads.used, curr.reads.cluster, curr.core, curr.low, curr.clust, curr.pre.clust,
               curr.dropped.lack, curr.dropped.low, curr.dropped.fail, curr.dropped.low.clone,
               curr.pcr.correct, curr.clone.dropped.lowq, curr.clone.preclust)
    assembleNames_v <- names(assembleRow_v)
  
  # Populate data frame
  align_df <- rbind(align_df, alignRow_v, stringsAsFactors = F)
  colnames(align_df) <- alignNames_v
  assemble_df <- rbind(assemble_df, assembleRow_v, stringsAsFactors = F)
  colnames(assemble_df) <- assembleNames_v

}  #  for

### Output
write.table(align_df, 
            file=file.path(outDir_v, "mixcr.alignment.QC.summary.txt"),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

write.table(assemble_df,
	    file=file.path(outDir_v, "mixcr.assembly.QC.summary.txt"),
	    quote=FALSE,
	    sep="\t",
	    row.names=FALSE)
