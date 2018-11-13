#!/usr/bin/Rscript

###
### Subset clones by a specific frequency group or groups
###

### Take normalized clone count files and subset them based on the information in the freqGroups file

### Dependencies
suppressMessages(library(data.table))
suppressMessages(library(optparse))
source("/home/exacloud/lustre1/users/hortowe/2016_11_27_stable_repos/WesPersonal/utilityFxns.R")

### TODO - should change outFile to outDir and write out the complete file, but also each frequency division as well

### Make list of options
optlist <- list(
  make_option(
    c("-i", "--normDir"),
    type = "character",
    help = "Directory of MiXCR clone files. Should be collapsed!"),
  make_option(
    c("-f", "--freqFile"),
    type = "character",
    help = "File of clones subset by frequency group. Preferably the 'full' file, but can be a subset if --group fits."),
  make_option(
    c("-d", "--div"),
    type = "character",
    help = "comma-separated, no space list of frequency groups to subset. e.g. 'Hyperexpanded,Large,Rare'"),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    help = "Directory to output files"),
  make_option(
    c("-l", "--log"),
    type = "logical",
    help = "TRUE if log file should be written to outDir. FALSE if no log should be written.")
)

### Parse commandline
p <- OptionParser(usage = "%prog -i normDirectory -f freqFile -o outputDirectory -l log(T/F)",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

normDir_v <- args$normDir
freqFile_v <- args$freqFile
groups_v <- args$div
outDir_v <- args$outDir
log_v <- args$log

### Print log
if (log_v){
    returnSessionInfo(args_lsv = args, out_dir_v = outDir_v)
} # fi

### Handle group argument
groups_v <- unlist(strsplit(groups_v, split = ","))

### Get files and names
normFiles_v <- list.files(normDir_v)
normFiles_v <- normFiles_v[order(as.numeric(gsub("^.*_S|_.*", "", normFiles_v)))]
batchName_v <- strsplit(normFiles_v[1], split = "_")[[1]][1]

### Read in freq-group data
freqData_dt <- fread(freqFile_v)
#freqName_v <- gsub(batchName_v, "", basename(freqFile_v)); freqName_v <- gsub("_|_clones.txt", "", freqName_v)

### Get freq column names
freqCols_v <- grep("nb|Normalized", colnames(freqData_dt), value = T)

### Set merge columns
mergeCols_v <- c("clonalSequence", "aaSeqCDR3", "V segments", "J segments", freqCols_v)
#mergeCols_v <- c("Clonal sequence(s)", "AA. Seq. CDR3", "V segments", "J segments", freqCols_v)

### Set columns to remove
toRem_v <- c("cloneId", "index")
#toRem_v <- c("Clone ID")

### For each sample, grab the appropriate clones and then export
for (i in 1:length(normFiles_v)){

	## Get name of file and sample name
	currFile_v <- normFiles_v[i]
	currSample_v <- grep("S[0-9]+", unlist(strsplit(currFile_v, split = "_")), value = T)
	#print(currSample_v)

	## Get data and columns
	currData_dt <- fread(file.path(normDir_v, currFile_v))
	currColumns_v <- colnames(currData_dt)

	## Remove duplicate columns
	dupCols_v <- which(duplicated(currColumns_v))
	if (length(dupCols_v) > 0) {
		currData_dt <- currData_dt[,-dupCols_v, with = F]
		currColumns_v <- colnames(currData_dt)
	} # fi

	## Subset Freq data
	currFreqSub_dt <- freqData_dt[Sample == currSample_v & Div %in% groups_v,]

	## Get the same from norm data
	currMerge_dt <- merge(currData_dt, currFreqSub_dt, by = mergeCols_v, all = F)

	## Subset back to original columns
	currMerge_dt <- currMerge_dt[,mget(currColumns_v)]

	## Skip empty
	if (currMerge_dt[,.N] == 0){
		next
	} # fi

	## Remove ID columns
	currMerge_dt <- currMerge_dt[,(toRem_v) := NULL]

	## Fix semicolons from collapsed clones
	currMerge_mat <- apply(currMerge_dt, c(1,2), function(x) {
			## Not sure if element is dt or vector
			x1 <- unlist(x)
			## Remove edge case where semicolon at beginning
			x1 <- gsub("^;", "", x1)
			## Split and get unique
			x1 <- unique(strsplit(x1, split = ";")[[1]])
			## Handle edge case of blank
			if (length(x1) == 0) {x1 <- ''}
			## Trim WS
			x1 <- trimws(x1)
			## Turn into number, if needed
			if (suppressWarnings(is.na(unique(suppressWarnings(as.numeric(x1)))))) {
				x1 <- x1[1]
				
			} else {
				x1 <- as.numeric(x1)
				x1 <- sum(x1)
			}
			return(x1)})
	
	## Write output
	freqName_v <- paste0(groups_v, collapse = '')
	currOutName_v <- paste(paste(batchName_v, currSample_v, freqName_v, sep = "_"), "txt", sep = ".")
	write.table(currMerge_mat, file = file.path(outDir_v, currOutName_v), sep = '\t', quote = F, row.names = F)
} # for i
