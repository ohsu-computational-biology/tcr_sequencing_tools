#!/usr/bin/Rscript

###
### Count Singletons
###

### For each sample in a TCR-seq batch, find out how many have a clone count of 1.

####################
### DEPENDENCIES ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

library(wrh.rUtils)
library(optparse)

####################
### COMMAND LINE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

optlist <- list(
  make_option(
    c("-i", "--inputDir"),
    type = "character",
    help = "Directory of normalized MiXCR clone files"
  ),
  make_option(
    c("-c", "--column"),
    type = "character",
    default = "nb.clone.count",
    help = "Name of column to be used. Should be nb.clone.count (default) or cloneCount" 
),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    help = "Output directory"
  )
)

p <- OptionParser(usage = "%prog -i inputDir -g groups -r rest -o outDir",
                  option_list = optlist)

args <- parse_args(p)
opt <- args$options

inputDir_v <- args$inputDir
column_v <- args$column
outDir_v <- args$outDir

############
### BODY ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############

### Read in
inputFiles_v <- list.files(inputDir_v)
inputFiles_v <- inputFiles_v[order(as.numeric(gsub("^.*_S|_.*|\\..*$", "", inputFiles_v)))]

### Empty vars
singletons_v <- names_v <- NULL

### Run for each
for (i in 1:length(inputFiles_v)) {
  
  ## Get name and read in
  currFile_v <- inputFiles_v[i]
  currName_v <- paste0("S", gsub("^.*_S|\\..*$|_align.*", "", currFile_v))
  currData_dt <- fread(file.path(inputDir_v, currFile_v))
  
  ## Count number of singletons
  currN_v <- nrow(currData_dt[get(column_v) == 1,])
  
  ## Add
  singletons_v <- c(singletons_v, currN_v)
  names_v <- c(names_v, currName_v)
  
} # for i

### Output
out_dt <- data.table("Sample" = names_v, "nSingleton" = singletons_v)

### Get batch
batch_v <- unique(gsub("_S[0-9]+.*", "", inputFiles_v))
batch_v <- gsub("_.*", "", batch_v)

### Write output
write.table(out_dt,
            file = file.path(outDir_v, paste0(batch_v, "_singletons.txt")),
            sep = "\t", quote = F, row.names = F)
