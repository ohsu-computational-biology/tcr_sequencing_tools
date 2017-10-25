#!/usr/bin/Rscript

###
### Collapse clones with identical V, J, and AA CDR3, but different NT CDR3. Sum counts.
###

### Dependencies
library(data.table)
library(optparse)
source("/home/exacloud/lustre1/users/hortowe/2016_11_27_stable_repos/WesPersonal/utilityFxns.R")

### TODO - should change outFile to outDir and write out the complete file, but also each frequency division as well

### Make list of options
optlist <- list(
  make_option(
    c("-i", "--inputDir"),
    type = "character",
    help = "Directory of MiXCR clone files"),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    help = "Directory to output new clone files"),
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

inputDir_v <- args$inputDir
outDir_v <- args$outDir
qcDir_v <- args$qcDir
log_v <- args$log

### Print log
if (log_v){
    returnSessionInfo(args_lsv = args, out_dir_v = outDir_v)
} # fi

### Get files and names
inputFiles_v <- list.files(inputDir_v)
inputFiles_v <- inputFiles_v[order(as.numeric(gsub("^.*_S|_.*", "", inputFiles_v)))]
inputNames_v <- sapply(inputFiles_v, function(x) grep("S[0-9]+", unlist(strsplit(x, split = "_")), value = T), USE.NAMES = F)
batchName_v <- strsplit(inputFiles_v[1], split = "_")[[1]][1]

### Read in data
clones_lsdt <- sapply(inputFiles_v, function(x) {
    ## Get data
    y <- fread(file.path(inputDir_v, x), drop = "reads")
    return(y)}, simplify = F)

names(clones_lsdt) <- inputNames_v

### Remove any empty files
clones_lsdt <- clones_lsdt[sapply(clones_lsdt, function(x) dim(x)[1]) > 0]

### Collect number of clones
origCloneCount_v <- sapply(clones_lsdt, nrow)

### Combine Clones

seqCol_v <- grep("aa.*CDR3|AA.*CDR3", colnames(clones_lsdt[[1]]), value = T)
groupCols_v <- c("V segments", seqCol_v, "J segments")
aggCols_v <-  colnames(clones_lsdt[[1]])[!(colnames(clones_lsdt[[1]]) %in% groupCols_v)]

clones_lsdt <- sapply(clones_lsdt, function(x) {
    ## Collapse like clones
    new_dt <- x[, lapply(.SD, function(y) paste(y, collapse = ';')), by = groupCols_v]
    ## Add back clone counts and freqs
    new_dt$`Normalized clone count` <- sapply(new_dt$`Normalized clone count`, function(y) 
						sum(as.numeric(unlist(strsplit(y, split = ";")))), USE.NAMES = F)
    
    new_dt$`nb.clone.count` <- sapply(new_dt$`nb.clone.count`, function(y) 
						sum(as.numeric(unlist(strsplit(y, split = ";")))), USE.NAMES = F)
    
    new_dt$`Normalized clone fraction` <- sapply(new_dt$`Normalized clone fraction`, function(y) 
						sum(as.numeric(unlist(strsplit(y, split = ";")))), USE.NAMES = F)
    
    new_dt$`nb.clone.fraction` <- sapply(new_dt$`nb.clone.fraction`, function(y) 
						sum(as.numeric(unlist(strsplit(y, split = ";")))), USE.NAMES = F)
    ## Check fractions
    if (!(all.equal(sum(new_dt$`Normalized clone fraction`), 1))) stop("Incorrect 'Normalized clone fraction' summation")
    if (!(all.equal(sum(new_dt$`nb.clone.fraction`), 1))) stop("Incorrect 'nb.clone.fraction' summation")
    ## Return
    return(new_dt)}, simplify = F)

### Collect new number of clones
newCloneCount_v <- sapply(clones_lsdt, nrow)
    
### Combine new and old
cloneCountSummary_mat <- rbind(origCloneCount_v, newCloneCount_v)
cloneCountSummary_dt <- as.data.table(t(cloneCountSummary_mat))
colnames(cloneCountSummary_dt) <- c("original", "collapsed")

### Prepare output names
outputFiles_v <- sapply(inputFiles_v, function(x) gsub("\\.txt", "_collapsed\\.txt", x), USE.NAMES = F)

### Write out tables
for (i in 1:length(outputFiles_v)){
    ## Get file
    currData_dt <- clones_lsdt[[i]]
    ## Get name
    currInNum_v <- names(clones_lsdt)[i]
    ## Get other name
    currOutNum_v <- grep("S[0-9]+", unlist(strsplit(outputFiles_v[i], split = "_")), value = T)
    ## Write
    if (all.equal(currInNum_v, currOutNum_v)) {
        write.table(currData_dt,
            file.path(outDir_v, outputFiles_v[i]), sep = '\t', quote = F, row.names = F)
    } # fi
} # for

### Write QC
write.table(cloneCountSummary_dt, file = file.path(qcDir_v, "cloneCollapseSummary.txt"), sep = '\t', quote = F, row.names = F)

