#!/usr/bin/Rscript

###
### Top Clone Frequencies
###

### For each sample in a TCR-seq batch, sort by clone frequency and then find the sum of
### the top N clones' frequencies.

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
    c("-g", "--groups"),
    type = "character",
    help = "List of top clones to take sums of. Comma-separated no spaces. (e.g. 5,10,50,100). Can also specify ranges (e.g. 6-25)"
    ),
  make_option(
    c("-r", "--rest"),
    type = "logical",
    help = "logical. TRUE - include all remaining clones as a group, e.g. 101-the rest. FALSE - stop at largest group above"
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
groups_v <- splitComma(args$groups)
rest_v <- args$rest
outDir_v <- args$outDir

############
### BODY ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############

### Read in
inputFiles_v <- list.files(inputDir_v)
inputFiles_v <- inputFiles_v[order(as.numeric(gsub("^.*_S|_.*|\\..*$", "", inputFiles_v)))]

### Add rest
if (rest_v) groups_v <- c(groups_v, "rest")

### Make empty matrix
out_mat <- matrix(nrow = length(inputFiles_v), ncol = length(groups_v))
colnames(out_mat) <- groups_v
rownames(out_mat) <- 1:length(inputFiles_v)

### Run for each
for (i in 1:length(inputFiles_v)) {
  
  ## Get name and read in
  currFile_v <- inputFiles_v[i]
  currName_v <- paste0("S", gsub("^.*_S|\\..*$|_align.*", "", currFile_v))
  currData_dt <- fread(file.path(inputDir_v, currFile_v))
  
  ## Sort
  currData_dt <- currData_dt[order(nb.clone.fraction, decreasing = T)]
  
  ## Sum each freq
  for (j in 1:length(groups_v)) {
    
    ## Get group
    currGroup_v <- groups_v[j]

    ## Special if 'rest'
    if (currGroup_v == 'rest') {
	one_v <- as.numeric(gsub("^.*\\-", "", groups_v[j-1]))+1
	two_v <- nrow(currData_dt)
	sum_v <- sum(currData_dt[one_v:two_v, nb.clone.fraction])
    } else if (length(grep("\\-", currGroup_v)) > 0) {
	one_v <- as.numeric(gsub("\\-.*$", "", currGroup_v))
	two_v <- as.numeric(gsub("^.*\\-", "", currGroup_v))
	two_v <- ifelse(two_v <= nrow(currData_dt), two_v, nrow(currData_dt))
	sum_v <- sum(currData_dt[one_v:two_v, nb.clone.fraction])
    } else {
	val_v <- as.numeric(currGroup_v)
	val_v <- ifelse(val_v <= nrow(currData_dt), val_v, nrow(currData_dt))
        sum_v <- sum(currData_dt[1:val_v, nb.clone.fraction])
    } # fi

    ## Add
    out_mat[i,currGroup_v] <- sum_v
    
  } # for j
  
  rownames(out_mat)[i] <- currName_v
  
} # for i

### Format output
out_dt <- convertDFT(out_mat, newName_v = "Sample")

### Get batch
batch_v <- unique(gsub("_S[0-9]+.*", "", inputFiles_v))

### Write output
write.table(out_dt,
            file = file.path(outDir_v, paste0(batch_v, "_topCloneFreq.txt")),
            sep = "\t", quote = F, row.names = F)
