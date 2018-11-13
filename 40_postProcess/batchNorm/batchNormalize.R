#!/usr/bin/Rscript

###
### Batch to Batch Normalization prodecure
###

### Perform median-normalization across samples usinge the formula:
  ### Adjusted(ij) = Raw(ij) * median[C(1k)...C(jk)]/Cj
  ### Where i is an individual clone in a sample
  ### j is a sample
  ### k is a sequencing batch
  ### Cj is the sum of all clone counts in sample j
  ### medain[C(1k)...C(jk)] is the median of the sums of all samples 1..j in all batches k

### Dependencies
suppressMessages(library(data.table))
library(optparse)
suppressMessages(library(affy))
suppressMessages(library(mcr))
library(ggplot2)
source("/home/exacloud/lustre1/CompBio/users/hortowe/2016_11_27_stable_repos/WesPersonal/utilityFxns.R")
source("/home/exacloud/lustre1/CompBio/users/hortowe/tcr_sequencing_tools/normalize/batchNormHelperFxn.R")


### Make list of options
optlist <- list(
  make_option(
    c("-b", "--baseDir"),
    type = "character",
    help = "Base directory that contains all pertinent directories and data for run."
  ),
  make_option(
    c("-d", "--dataDir"),
    type = "character",
    help = "sub-path to clone files within the baseDir. This directory contains all samples from all batches."
  ),
  make_option(
    c("-m", "--meta"),
    type = "character",
    help = "Metadata file containing treatment designations at a minimum. 
    Column1 = Sample; Column2 = Treatment; Column3 = Tissue; Column4 = Batch"
  ),
  make_option(
    c("-n", "--norm"),
    type = "character",
    help = "Character to determine which count column to use. 'raw' = raw counts; 'orig' = orig norm method; 'nb' = new, negative binomial method'."
  ),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    help = "sub-path to directory for writing output files, within baseDir."
  ),
  make_option(
    c("-p", "--pre"),
    type = "character",
    help = "prefix and suffix to add around batch name when outputting files. comma-sep, no spaces. If two different batches, separte them by '.'"
  ),
  make_option(
    c("-f", "--debug"),
    type = "logical",
    default = FALSE,
    help = "Logical. TRUE - print session info and extra output to help with debugging. Also do not write output (tables and images). FALSE - normal output and write output files (tables and images)."
  ),
  make_option(
    c("-l", "--log"),
    type = "logical",
    default = FALSE,
    help = "Logical. TRUE - output session info. FALSE - do not output session info. If debug is TRUE and log is FALSE, then session info will be printed to STDOUT. If neither are set, no output."
  )
)


### Parse commandline
p <- OptionParser(usage = "%prog -1 inputDir1 -2 inputDir2 -m metaDataFile -o outputDirectory -d debug -l log",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

### Commands
baseDir_v <- args$baseDir
dataDir_v <- args$dataDir
metaFile_v <- args$meta
outDir_v <- args$outDir
norm_v <- args$norm
debug_v <- args$debug
log_v <- args$log
fix_v <- args$pre


##################
### PREPROCESS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################

### Read in metadata
metadata_dt <- fread(metaFile_v)

### Get batches and treatments
batchCol_v <- grep("atch", colnames(metadata_dt), value = T)
treatCol_v <- grep("reat", colnames(metadata_dt), value = T)
tissueCol_v <- grep("issue", colnames(metadata_dt), value = T)
sampleCol_v <- grep("ample", colnames(metadata_dt), value = T)

batches_v <- unique(metadata_dt[[batchCol_v]])
treatments_v <- unique(metadata_dt[[treatCol_v]])

### Read in data
### List where each element is a list containing the clone data.tables for one batch
### Name of each base element is the batch name. The name of that list's elements are the sample numbers
clones_lslsdt <- readDir_clone(dataDir = baseDir_v, subDir = dataDir_v, meta = metadata_dt, norm_v = norm_v)
print("Done reading in data.")

### Clean Data
### Remove pseudo clones and adjust V121/V122
clones_lslsdt <- sapply(clones_lslsdt, function(x) sapply(x, function(y) cln_clones(y), simplify = F), simplify = F)
print("Done cleaning data.")

#####################
### NORMALIZATION ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####################

### Get sum of clone counts for each sample
batchSums_v <- unlist(lapply(clones_lslsdt, function(x) unlist(lapply(x, function(y) sum(y$CloneCount)))))

### Get median
batchMedian_v <- median(batchSums_v)

### Normalize samples
normClones_lslsdt <- lapply(clones_lslsdt, function(x) lapply(x, function(y) norm4lib(clonedf = y, median_v = batchMedian_v)))

### Add Frequency
normClones_lslsdt <- lapply(normClones_lslsdt, function(x) lapply(x, function(y) {
	sum_v <- sum(y$normC)
	freq_v <- y$normC / sum_v
	y$normFreq <- freq_v
	return(y)}))
print("Done normalizing data.")


##############
### OUTPUT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############

### Make sure directory exists
outDir_v <- mkdir(baseDir_v, outDir_v)

### Split prefix and suffix
fix_v <- strsplit(fix_v, split = ",")[[1]]

for (i in 1:length(normClones_lslsdt)) {
  currBatch_v <- names(normClones_lslsdt)[i]
  currData_lsdt <- normClones_lslsdt[[i]]
  for (j in 1:length(currData_lsdt)) {
    currSample_v <- names(currData_lsdt)[j]
    currData_dt <- currData_lsdt[[j]]
    currOutName_v <- paste0(fix_v[1], currBatch_v, fix_v[2], "_S", currSample_v, "_batchNorm_clones.txt")
    write.table(currData_dt, file.path(outDir_v, currOutName_v), sep = '\t', quote = F, row.names = F)
  } # for j
} # for i


############################
### TISSUE NORMALIZATION ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############################

# ### Trying a different method where I take medians for each tissue
# tissues_lsv <- sapply(batches_v, function(x) unique(metadata_dt[get(batchCol_v) == x, get(tissueCol_v)]))
# names(tissues_lsv) <- batches_v
# 
# ### Get medians
# tissueMedians_v <- sapply(names(tissues_lsv), function(x) {
#   sapply(tissues_lsv[[x]], function(y){
#     ## Get samples
#     samples_v <- as.character(metadata_dt[get(batchCol_v) == x &
#                                             get(tissueCol_v) == y, get(sampleCol_v)])
#     samples2_v <- metadata_dt[get(batchCol_v) == x &
#                                 get(tissueCol_v) == y, get(sampleCol_v)]
#     ## Get sums
#     sums_v <- sapply(clones_lslsdt[[x]][samples_v], function(z) {
#       sum(z$CloneCount)
#     })
#     ## Get median
#     med_v <- median(sums_v)
#   })
# })
# 
# ### Get median of medians
# tissueMedian_v <- median(unlist(tissueMedians_v))
# 
# ### Normalize samples
# tissueNorm_lslsdt <- lapply(clones_lslsdt, function(x) lapply(x, function(y) norm4lib(clonedf = y, median_v = tissueMedian_v)))
