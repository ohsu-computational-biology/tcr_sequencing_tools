#!/usr/bin/Rscript
### RScript that aggregates reports from MiXCR's alignment tool
### Basically turn their report files into a single text file that can be easily manipulated

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
inputFiles_v <- inputFiles_v[order(as.numeric(gsub("^S|_align.*", "", inputFiles_v)))]

### Create empty data.frame
output_df <- data.frame()

### Extract info from each report
for (i in 1:length(inputFiles_v)){

    ## Get Record
    currFile_v <- inputFiles_v[i]
    curr.record <- readLines(file.path(inputDir_v, currFile_v))

    ## Update user
    if (length(curr.record) != 17){
	print(sprintf("Sample %s has %d lines. The norm is 17.", currFile_v, length(curr.record)))
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
    curr.Alignment.Aided <- extractRecord("Alignment-aided", curr.record, 
					c("num.align.aided.overlap", "pct.align.aided.overlap"))
    curr.Overlap.Not.Align <- extractRecord("Overlapped and not", curr.record,
                                            c("num.overlapped.and.not.algned", "pct.overlapped.and.not.aligned"))
    curr.TRB.chains <- extractRecord("TRB chains", curr.record, c("num.TRB.chains", "pct.TRB.chains"))

  ## Combine into a row to add to data frame
    row_v <- c(curr.date, curr.input, curr.output, curr.version, curr.total,
               curr.Success, curr.NoHit, curr.NoJ, curr.Overlap, curr.Overlap.Align, curr.Overlap.Not.Align)
    colNames_v <- names(row_v)
  
  # Populate data frame
  output_df <- rbind(output_df, row_v, stringsAsFactors = F)
  colnames(output_df) <- colNames_v

}  #  for

### Output
write.table(output_df, 
            file=file.path(outDir_v, "mixcr.alignment.QC.summary.txt"),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

