#!/usr/bin/Rscript

###
### Collapse clones with identical V, J, and AA CDR3, but different NT CDR3. Sum counts.
###

### Dependencies
library(data.table)
library(optparse)
source("/home/exacloud/lustre1/CompBio/users/hortowe/2016_11_27_stable_repos/WesPersonal/utilityFxns.R")

### TODO - should change outFile to outDir and write out the complete file, but also each frequency division as well

### Make list of options
optlist <- list(
  make_option(
    c("-i", "--inputFile"),
    type = "character",
    help = "Path to single MiXCR clone file"),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    help = "Path to and name of new clone file"),
  make_option(
    c("-q", "--qcDir"),
    type = "character",
    help = "Directory to output record of clones collapsed"),
  make_option(
    c("-l", "--log"),
    type = "logical",
    default = F,
    help = "TRUE if log file should be written to outDir. FALSE if no log should be written.")
)

### Parse commandline
p <- OptionParser(usage = "%prog -i inputDirectory -o outputDirectory -c columns -m metadata -d freqDivisions -l T/F",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

inputFile_v <- args$inputFile
outDir_v <- args$outDir
qcDir_v <- args$qcDir
log_v <- args$log

### Print log
if (log_v){
    returnSessionInfo(args_lsv = args, out_dir_v = outDir_v)
} # fi

### Read in data
#inputData_dt <- fread(inputFile_v, drop = "reads") # new format
inputData_dt <- fread(inputFile_v, drop = "Reads") # old format

sampleName_v <- grep("S[0-9]+", unlist(strsplit(inputFile_v, split = "_")), value = T)

### Collect number of clones
origCloneCount_v <- nrow(inputData_dt)

### Get columns
seqCol_v <- grep("aa.*CDR3|AA.*CDR3", colnames(inputData_dt), value = T)
groupCols_v <- c("V segments", seqCol_v, "J segments")
aggCols_v <-  colnames(inputData_dt)[!(colnames(inputData_dt) %in% groupCols_v)]

### Combine clones
new_dt <- inputData_dt[, lapply(.SD, function(y) paste(y, collapse = ';')), by = groupCols_v]
    
### Add back clone counts and freqs
new_dt$`Normalized clone count` <- sapply(new_dt$`Normalized clone count`, function(y) 
						sum(as.numeric(unlist(strsplit(y, split = ";")))), USE.NAMES = F)
    
new_dt$`Normalized clone fraction` <- sapply(new_dt$`Normalized clone fraction`, function(y) 
						sum(as.numeric(unlist(strsplit(y, split = ";")))), USE.NAMES = F)
### Check fractions
if (!(all.equal(sum(new_dt$`Normalized clone fraction`), 1))) stop("Incorrect 'Normalized clone fraction' summation")


### Do for nb, if present
if (length(grep("nb", colnames(new_dt))) > 0){
    new_dt$`nb.clone.fraction` <- sapply(new_dt$`nb.clone.fraction`, function(y) 
						sum(as.numeric(unlist(strsplit(y, split = ";")))), USE.NAMES = F)

    new_dt$`nb.clone.count` <- sapply(new_dt$`nb.clone.count`, function(y) 
						sum(as.numeric(unlist(strsplit(y, split = ";")))), USE.NAMES = F)
    
    if (!(all.equal(sum(new_dt$`nb.clone.fraction`), 1))) stop("Incorrect 'nb.clone.fraction' summation")
} # fi

### Collect new number of clones
newCloneCount_v <- nrow(new_dt)
    
### Combine new and old
cloneCountSummary_dt <- data.table("Sample" = sampleName_v, "original" = origCloneCount_v, "collapsed" = newCloneCount_v)

### Prepare output name
outputName_v <- gsub("\\.txt", "_collapsed\\.txt", basename(inputFile_v))

write.table(new_dt, file.path(outDir_v, outputName_v), sep = '\t', quote = F, row.names = F)

### Write QC
qcFile_v <- file.path(qcDir_v, "cloneCollapseSummary.txt")

if (file.exists(qcFile_v)) {append_v <- T; colNames_v <- F} else {append_v <- F; colNames_v = T}

write.table(cloneCountSummary_dt, file = qcFile_v, sep = '\t', quote = F, row.names = F, append = append_v, col.names = colNames_v)

