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
    c("-f", "--freqs"),
    type = "character",
    help = "List of top clones to take sums of. Comma-separated no spaces. (e.g. 5,10,50,100)"
    ),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    help = "Output directory"
  )
)

p <- OptionParser(usage = "%prog -i inputDir -f freqs -o outDir",
                  option_list = optlist)

args <- parse_args(p)
opt <- args$options

inputDir_v <- args$inputDir
freqs_v <- splitComma(args$freqs)
outDir_v <- args$outDir

inputDir_v <- "/Users/hortowe/OHSU/tcr_spike/data/LIB190701LC/normClones"
freqs_v <- splitComma("5,10,15,20,50,100,200")

############
### BODY ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############

### Read in
inputFiles_v <- list.files(inputDir_v)
inputFiles_v <- inputFiles_v[order(as.numeric(gsub("^.*_S|_.*|\\..*$", "", inputFiles_v)))]

### Make empty matrix
out_mat <- matrix(nrow = length(inputFiles_v), ncol = length(freqs_v))
colnames(out_mat) <- freqs_v
rownames(out_mat) <- 1:length(inputFiles_v)

### Run for each
for (i in 1:length(inputFiles_v)) {
  
  ## Get name and read in
  currFile_v <- inputFiles_v[i]
  currName_v <- paste0("S", gsub("^.*_S|\\..*$", "", currFile_v))
  currData_dt <- fread(file.path(inputDir_v, currFile_v))
  
  ## Sort
  currData_dt <- currData_dt[order(nb.clone.fraction, decreasing = T)]
  
  ## Sum each freq
  for (j in 1:length(freqs_v)) {
    
    currFreq_v <- freqs_v[j]
    sum_v <- sum(currData_dt[1:currFreq_v, nb.clone.fraction])
    out_mat[i,currFreq_v] <- sum_v
    
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
