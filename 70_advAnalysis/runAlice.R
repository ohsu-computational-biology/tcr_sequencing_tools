#!/usr/bin/Rscript

###
### RUN ALICE
###

####################
### DEPENDENCIES ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

library(wrh.rUtils)
library(optparse)

####################
### COMMAND LINE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

### Make list of options
optlist <- list(
  make_option(
    c("-i", "--inputDir"),
    type = "character",
    help = "Directory containing TCR files in ALICE format."),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    help = "path to directory to write output files"),
  make_option(
    c("-l", "--loadDir"),
    type = 'character',
    default = "/home/exacloud/lustre1/CompBio/users/hortowe/myApps/ALICE",
    help = "path to ALICE github repo location."),
  make_option(
    c("-c", "--cpu"),
    type = "numeric",
    default = 1,
    help = "Number of cpus to run on"),
  make_option(
    c("-n", "--nSeq"),
    type = "numeric",
    default = 1000000,
    help = "Number of sequences...not sure where this is applied"
  ),
  make_option(
    c("-t", "--treat"),
    type = "logical",
    default = F,
    help = "TRUE - one file for each treatment; FALSE - one directory for each treatment"
  )
)

### Parse commandline
p <- OptionParser(usage = "%prog -i inputDir -o outDir -l loadDir -c cpu -n nSeq -t treat",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

### Assign
inputDir_v <- args$inputDir
outDir_v <- args$outDir
loadDir_v <- args$loadDir
cpu_v <- args$cpu
n_v <- args$nSeq
byTreat_v <- args$treat

print(inputDir_v)
print(outDir_v)
print(loadDir_v)

############
### BODY ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############

### Set directory and load script
setwd(loadDir_v)
source("ALICE.R")

### Set output
if (!dir.exists(outDir_v)) mkdir(dirname(outDir_v), basename(outDir_v))
setwd(outDir_v)

### Get files
inputFiles_v <- list.files(inputDir_v)
index_v <- ifelse(byTreat_v, 1, 2)
inputNames_v <- sapply(inputFiles_v, function(x) strsplit(x, split = "_")[[1]][index_v])

### Read in
inputData_lsdt <- lapply(inputFiles_v, function(x) fread(file.path(inputDir_v, x)))
names(inputData_lsdt) <- inputNames_v

### Evaluate - this outputs some results
outData_lsdt <- ALICE_pipeline(DTlist                 = inputData_lsdt,
                               folder                 = "rda",
                               cores                  = cpu_v,
                               iter                   = 10,
                               nrec                   = n_v,
			       P_thres                = 0.001,
			       cor_method             = "BH",
			       qL                     = F,
			       Read_count_filter      = 0,
			       Read_count_neighbour    = 1)

### Notify user of any empties
for (i in 1:length(outData_lsdt)) {
	currName_v <- names(outData_lsdt)[i]
	currData_dt <- outData_lsdt[[currName_v]]
	if (nrow(currData_dt) == 0) {
		cat(sprintf("No significant hits found for: %s\n", currName_v))
	} # fi
} # for

### Write output
outDir_v <- mkdir(outDir_v, "results")
lapply(names(outData_lsdt), function(x) write.table(outData_lsdt[[x]], file = file.path(outDir_v, paste0(x, ".txt")), sep = '\t', quote = F, row.names = F))
