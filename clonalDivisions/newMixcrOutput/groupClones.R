#!/usr/bin/Rscript

###
### Group Clones by different frequency classifications
###

### Take a batch of clone files, combine all of the files togther, and group them by their different frequencies

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
    help = "output file containing all samples plus frequency designations"),
  make_option(
      c("-c", "--columns"),
      type = "character",
      default = c("Normalized clone count,Normalized clone fraction,clonalSequence,aaSeqCDR3,V segments,J segments"),
      help = "Column names to read in. Count and fraction are required, but all default columns are recommended. Comma-separated, no spaces.\
		Example: 'Col1,Col 2,Col 3 a'"),
    make_option(
    c("-m", "--meta"),
    type = "character",
    help = "Metadata file containing at a minimum treatment designations for each sample."),
  make_option(
      c("-d", "--divisions"),
      type = "character",
      default = c("Rare = 0.00001,Small = 0.0001,Medium = 0.001,Large = 0.01,Hyperexpanded = 1"),
      help = "Clonal frequency divisions to group clones into. Comma-separated, no spaces. Requires name and value separated by equals.\
		Example: 'Rare = 0.00001,Small=0.0001,Medium = 0.001'"),
    make_option(
    c("-l", "--log"),
    type = "logical",
    help = "TRUE if log file should be written to outDir. FALSE if no log should be written.")
)

### Parse commandline
p <- OptionParser(usage = "%prog -i inputDirectory -o outputDirectory -c columns -m metadata -d freqDivisions -l T/F",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

inputDir_v <- args$inputDir
metaFile_v <- args$meta
outDir_v <- args$outDir
log_v <- args$log

### Handle non-traditional arguments
columns_v <- args$columns; columns_v <- unlist(strsplit(columns_v, split = ','))
divisions_v <- args$divisions; divisions_v <- unlist(strsplit(divisions_v, split = ','))
divNames_v <- sapply(divisions_v, function(x) trimws(strsplit(x, split = '=')[[1]][1]), USE.NAMES = F)
divNums_v <- sapply(divisions_v, function(x) trimws(strsplit(x, split = '=')[[1]][2]), USE.NAMES = F)
divisions_v <- as.numeric(divNums_v); names(divisions_v) <- divNames_v

### Print log
if (log_v){
    returnSessionInfo(args_lsv = args, out_dir_v = outDir_v)
} # fi

### Get files and names
inputFiles_v <- list.files(inputDir_v)
inputFiles_v <- inputFiles_v[order(as.numeric(gsub("^.*_S|_.*", "", inputFiles_v)))]
batchName_v <- strsplit(inputFiles_v[1], split = "_")[[1]][1]

### Get metadata
meta_dt <- fread(metaFile_v)
sampleCol_v <- grep("ample", colnames(meta_dt), value = T)
treatCol_v <- grep("eatment", colnames(meta_dt), value = T)

### Subset to only contain samples in meta
baseFile_v <- strsplit(inputFiles_v[1], split = "S[0-9]+")[[1]]
toKeep_v <- paste0(baseFile_v[1], "S", unlist(meta_dt[,get(sampleCol_v)]), baseFile_v[2])
inputFiles_v <- inputFiles_v[inputFiles_v %in% toKeep_v]
#inputNames_v <- sapply(inputFiles_v, function(x) strsplit(x, split = "_")[[1]][2], USE.NAMES = F)
inputNames_v <- sapply(inputFiles_v, function(x) grep("S[0-9]+", unlist(strsplit(x, split = "_")), value = T), USE.NAMES=F)

### Read in data
clones_lsdt <- sapply(inputFiles_v, function(x) {
    ## Get data
    y <- fread(file.path(inputDir_v, x), select = columns_v)
    ## Get sample character
    sample_v <- grep("S[0-9]+", unlist(strsplit(x, split = "_")), value = T)
    #sample_v <- strsplit(x, split = "_")[[1]][2]
    ## Get sample number
    sampNum_v <- as.numeric(gsub("S", "", sample_v))
    ## Add sample character to data
    y$Sample <- sample_v
    ## Add treatment to data
    y$Treatment <- meta_dt[get(sampleCol_v) == sampNum_v, get(treatCol_v)]
    ## Add clone index column
    y$id <- 1:nrow(y)
    ## Add empty division column
    y$Div <- character()
    return(y)}, simplify = F)

names(clones_lsdt) <- inputNames_v

print("Finished reading data")

### Remove any empty files
clones_lsdt <- clones_lsdt[sapply(clones_lsdt, function(x) dim(x)[1]) > 0]

### Classify Clones

clones_lsdt <- sapply(clones_lsdt, function(x) {
    ## Create empty column
#    x$Div <- character()
print(head(x))
    ## Iterate for each row and determine which divisions are greater than the current row's frequency
    ## Since the division values represent the upper limit of that grouping, the first division that is
    ## greater than the frequency is the one that we want. (e.g. Large is 0.001 to 0.01, Medium is 0.0001 to 0.001,
    ## and Hyper is 0.01 to 1. Given freq of 0.00099, which is less than all 3, medium will be chosen. Given freq of
    ## freq of 0.0010001, which is less than Large and Hyper, Large will be chosen. Given freq of 0.01, Large and Hyper
    ## will match, and large will be chosen.)
    x[, Div := sapply(x[,`Normalized clone fraction`], function(y) names(which(divisions_v >= y)[1]))]
print(head(x))
    return(x)
}, simplify = F)

### Combine into giant data.table
clones_dt <- do.call(rbind, clones_lsdt)

### Gather summary info
summary_lsdt <- sapply(names(divisions_v), function(x) clones_dt[Div == x, .N, by = Sample], simplify = F)
summary_dt <- mergeDTs(summary_lsdt, mergeCol_v = "Sample", keepCol_v = "N")

### Split into multiple data.tables
divisions_lsdt <- list("summary" = summary_dt, "full" = clones_dt)
for (i in 1:length(divisions_v)){
    ## Get divisions
    currDiv_v <- names(divisions_v)[i]
    ## Subset
    currSubset_dt <- clones_dt[Div == currDiv_v,]
    ## Add to list
    divisions_lsdt[[currDiv_v]] <- currSubset_dt
} # for i

### Write out all of the data.tables
for (i in 1:length(divisions_lsdt)){
    currOut_dt <- divisions_lsdt[[i]]
    currName_v <- names(divisions_lsdt)[i]
    write.table(currOut_dt,
                file.path(outDir_v, paste(batchName_v, currName_v, "clones.txt", sep = "_")),
                sep = '\t', quote = F, row.names = F)
} # for i


warnings()

