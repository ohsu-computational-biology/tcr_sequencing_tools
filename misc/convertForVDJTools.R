#!/usr/bin/Rscript

###
### Convert subsetByGroup.R output so that it can be used with VDJtools
###

suppressMessages(library(data.table))
library(optparse)

### Make list of options
optlist <- list(
  make_option(
    c("-i", "--inputDir"),
    type = "character",
    help = "Directory containing clone files output by subsetByGroup.R"
  ),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    help = "File path to directory for writing output files"
  )
)


### Parse commandline
p <- OptionParser(usage = "%prog -i inputDir -o outputDirectory",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

### Commands
inputDir_v <- args$inputDir
outDir_v <- args$outDir

inputFiles_v <- list.files(inputDir_v)

# inputDir_v <- "~/Desktop/OHSU/tcr_spike/data/LIB170920LC/collapsed/newNorm/testClones/"
# inputFiles_v <- "LIB170920LC_S1_Hyperexpanded.txt"
# compareFile_v <- "~/Desktop/OHSU/tcr_spike/data/LIB170920LC/collapsed/newNorm/testClones/1.txt"
# compareData_dt <- fread(compareFile_v)
# compareEntry_dt <- currData_dt[nSeqCDR3 == "TGTGCCAGCAGCCACGGGACTGGGGTCTATGAGCAGTACTTC",]

for (i in 1:length(inputFiles_v)){
  ## Get data and info
  currFile_v <- inputFiles_v[i]
  currSplit_v <- strsplit(currFile_v, split = "_|\\.")[[1]]
  currBatch_v <- currSplit_v[1]
  currType_v <- currSplit_v[3]
  currNum_v <- gsub("S", "", currSplit_v[2])
  outName_v <- paste0(currNum_v, ".txt")
  currData_dt <- fread(file.path(inputDir_v, currFile_v))
  print(sprintf("Currently on sample: %s", currFile_v))
  ## Subset 
  currSub_dt <- currData_dt[,mget(c("nb.clone.count", "nb.clone.fraction", "nSeqCDR3", "aaSeqCDR3", "V segments", "bestDHit", "J segments", "bestVAlignment", "bestDAlignment", "bestJAlignment"))]
  ## Fix V, D, J
  currSub_dt$bestDHit <- gsub("\\*00", "", currSub_dt$bestDHit)
  currSub_dt$`V segments` <- paste0("TRB", currSub_dt$`V segments`)
  currSub_dt$`J segments` <- paste0("TRB", currSub_dt$`J segments`)
  ## Extract Alignment stuff
  currSub_dt[, bestVAlignment := sapply(bestVAlignment, function(x) (as.numeric(strsplit(x, split = "\\|")[[1]][5]) - 1), USE.NAMES = F)]
  currSub_dt[, bestDStart := sapply(bestDAlignment, function(x) (as.numeric(strsplit(x, split = "\\|")[[1]][4])), USE.NAMES = F)]
  currSub_dt[, bestDAlignment := sapply(bestDAlignment, function(x) (as.numeric(strsplit(x, split = "\\|")[[1]][5]) - 1), USE.NAMES = F)]
  currSub_dt[, bestJAlignment := sapply(bestJAlignment, function(x) (as.numeric(strsplit(x, split = "\\|")[[1]][4])), USE.NAMES = F)]
  ## Change NA to -1 and blanks to '.'
  for (j in seq_len(ncol(currSub_dt))) set(currSub_dt, which(is.na(currSub_dt[[j]])),j,-1)
  for (j in seq_len(ncol(currSub_dt))) set(currSub_dt, which(currSub_dt[[j]] == ""),j,".")
  ## Reorder
  currSub_dt <- currSub_dt[,mget(c("nb.clone.count", "nb.clone.fraction", "nSeqCDR3", "aaSeqCDR3", "V segments", "bestDHit", "J segments", "bestVAlignment", "bestDStart", "bestDAlignment", "bestJAlignment"))]
  ## Change Names
  colnames(currSub_dt) <- c("count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j", "VEnd", "DStart", "DEnd", "JStart")
  ## Write output
  write.table(currSub_dt, file = file.path(outDir_v, outName_v), sep = '\t', quote = F, row.names = F)
} # for i
