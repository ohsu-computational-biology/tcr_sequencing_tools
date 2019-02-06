################################################
### SUBSET NORMALIZED CLONES BASED ON COUNTS ###
################################################

library(data.table)
library(optparse)

### Want to apply some sort of minimum count cut-off to our samples.

####################
### COMMAND LINE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

### Construct
optlist <- list(
  make_option(
    c("-i", "--inputDir"),
    type = "character",
    help = "Path to directory of normalized clone files to subset"
  ),
  make_option(
    c("-o", "--outputDir"),
    type = "character",
    help = "Path to directory to write output files"
  ),
  make_option(
    c("-c", "--cutOff"),
    type = "numeric",
    default = 1,
    help = "Value to use as cut-off. All clones with normalized count LESS THAN OR EQUAL TO cut-off will be removed."
  ),
  make_option(
    c("-q", "--qcDir"),
    type = "character",
    default = NULL,
    help = "Path to output directory for QC results. If NULL (default) will write to outputDir"
  ),
  make_option(
   c("-l", "--old") ,
   type = "logical",
   default = F,
   help = "TRUE - use old normalization columns; FALSE - use new normalization columns."
  )
)

### Parse
p <- OptionParser(usage = "%prog -i inputDir -o outputDir -c cutOff -q qcDir -l old",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

### Get args
inputDir_v <- args$inputDir
outputDir_v <- args$outputDir
cutOff_v <- args$cutOff
qcDir_v <- args$qcDir
old_v <- args$old

### For testing
inputDir_v <- "~/OHSU/tcr_spike/data/LIB180515LC/data/normalized_clones/"
cutOff_v <- 1

### Handle NULL
if (is.null(qcDir_v)) qcDir_v <- outputDir_v

###############
### WRANGLE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###############

### Read in files
inputFiles_v <- list.files(inputDir_v)

### Sort by sample number
inputFiles_v <- inputFiles_v[order(as.numeric(gsub("^.*_S|_.*$|\\..*$", "", inputFiles_v)))]

### Add names
inputNames_v <- paste0("S", gsub("^.*_S|_.*$|\\..*$", "", inputFiles_v))
names(inputFiles_v) <- inputNames_v

### Make output names
suffix_v <- paste0("_filter", cutOff_v, ".txt")
outputFiles_v <- gsub(".txt", suffix_v, inputFiles_v)

### Get batch name
batchName_v <- strsplit(inputFiles_v[1], split = "_")[[1]][1]

### Read in Data
inputData_lsdt <- sapply(inputFiles_v, function(x) fread(file.path(inputDir_v, x)), simplify = F)

### Determine count column to select
if (old_v) {
  count_v <- "Normalized clone count"
} else {
  count_v <- "nb.clone.count"
} # fi

### Get fraction columns to re-calculate
freqCols_v <- grep("[Ff]raction", colnames(inputData_lsdt[[1]]), value = T)
countCols_v <- grep("[Cc]ount", colnames(inputData_lsdt[[1]]), value = T)
cat(sprintf("Count - Frequency combinations that will be used. PLEASE DOUBLE CHECK THIS!!! \n\t%s\n",
            paste(paste(freqCols_v, countCols_v, sep = "_-_"), collapse = "\n\t")))

### Special case for raw data
if (!column_v %in% colnames(inputData_lsdt[[1]])){
  warning("Data are not normalized. Performing count filter on raw counts, which is not advisable.")
  column_v <- grep("cloneFraction|Clone fraction", colnames(inputData_lsdt[[1]]), value = T)
  count_v <- grep("cloneCount|Clone count", colnames(inputData_lsdt[[1]]), value = T)
} # fi

##############
### FILTER ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############

summary_mat <- matrix(nrow = length(inputData_lsdt), ncol = 4)

for (i in 1:length(inputData_lsdt)) {
  
  ## Get name and data
  currName_v <- names(inputData_lsdt)[i]
  currData_dt <- inputData_lsdt[[currName_v]]
  
  ## Record original rows
  origRows_v <- nrow(currData_dt)
  
  ## Apply cut-off and get new rows
  currData_dt <- currData_dt[get(count_v) > cutOff_v,]
  newRows_v <- nrow(currData_dt)
  
  ## Re-calculate frequencies
  for (j in 1:length(freqCols_v)) {
    
    ## Get current columns
    currCount_v <- countCols_v[j]
    currFreq_v <- freqCols_v[j]
    
    ## Get sum
    currSum_v <- sum(currData_dt[[currCount_v]])
    
    ## Re-calculate frequencies
    set(currData_dt, j = currFreq_v, value = currData_dt[[currCount_v]] / currSum_v)
    
    ## Check
    currCheckSum_v <- sum(currData_dt[[currFreq_v]])
    if (currCheckSum_v != 1) warning(sprintf("New frequency calculation for - %s - in sample - %s - does not sum to 1. Instead it is: %d\n", 
                                             currFreq_v, currName_v, currCheckSum_v))
  } # for j
  
  ## Add back to list
  inputData_lsdt[[currName_v]] <- currData_dt
  
  ## Calculate difference in number of clones
  currDiff_v <- origRows_v - newRows_v
  currDiffPct_v <- round(currDiff_v / origRows_v * 100, digits = 2)
  
  ## Summary output
  summary_mat[i,] <- c(origRows_v, newRows_v, currDiff_v, currDiffPct_v)
  
}

### Format summary
colnames(summary_mat) <- c("origRows", "filterRows", "numFiltered", "pctFiltered")
summary_df <- cbind("Sample" = inputNames_v, as.data.frame(summary_mat))

##############
### OUTPUT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############

### Write summary
summaryName_v <- file.path(qcDir_v, paste0(batchName_v, "_countFilter_", cutOff_v, "_summary.txt"))
write.table(summary_df, file = summaryName_v, sep = '\t', quote = F, row.names = F)

### Write clones
sapply(inputNames_v, function(x) {
  write.table(inputData_lsdt[[x]],
              file = file.path(outputDir_v, outputFiles_v[[x]]),
              sep = '\t', quote = F, row.names = F)
})
