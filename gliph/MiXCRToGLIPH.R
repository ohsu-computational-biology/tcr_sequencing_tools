#!/usr/bin/Rscript

###
### Convert MiXCR to GLIPH
###

### Convert clone files in the MiXCR format into GLIPH-compatible tables
### All samples from a batch will either be combined into one file or all samples from each treatment will be combined into one.

### Format:
  ### Column 1 : CDR3 AA seq
    ### MiXCR - "AA. Seq. CDR3"
    ### GLIPH - "CDR3b"
  ### Column 2 : V sequence
    ### MiXCR - "Best V hit"; format TRBV13-1*00
    ### GLIPH - "TRBV"; format TRBV13-1
  ### Column 3 : J sequence
    ### MiXCR - "Best J hit"; format TRBJ2-5*00
    ### GLIPH - "TRBJ"; format TRBJ2-5
  ### Column 4 : Patient (Sample, but want consistent names with tool)
    ### MiXCR - from file name
  ### Column 5 : Counts
    ### MiXCR - Normalized clone count


### Dependencies
library(data.table)
library(optparse)

### Make list of options
optlist <- list(
  make_option(
    c("-i", "--inputDir"),
    type = "character",
    help = "Directory of MiXCR clone files"
  ),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    help = "output directory to write GLIPH-compatible clone files"
  ),
  make_option(
    c("-m", "--meta"),
    type = "character",
    help = "If set, will output one file for each treatment, containing all samples in that treatment. If not, will output one file containing all samples in the batch.")
)

### Parse commandline
p <- OptionParser(usage = "%prog -i inputDirectory -o outputDirectory",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

### Commands
cloneDir_v <- args$inputDir
outDir_v <- args$outDir
metaFile_v <- args$meta

### Get files and meta
cloneFiles_v <- list.files(cloneDir_v)
if (!is.null(metaFile_v)){
    ## Read in
    meta_dt <- fread(metaFile_v)

    ## Get treatment and sample columns from metadata
    treatCol_v <- grep("eatment", colnames(meta_dt), value = T)
    sampleCol_v <- grep("ample", colnames(meta_dt), value = T)
} # fi   

### Get batch
batchName_v <- strsplit(cloneFiles_v[1], split = "_")[[1]][1]

### Column variables
toRead_v <- c("AA. Seq. CDR3", "Best V hit", "Best J hit", "Normalized clone count")
vj_v <- c("Best V hit", "Best J hit")

### Read in files and fix up
clones_lsdt <- sapply(cloneFiles_v, function(x) {
    ## Get data
    y <- fread(file.path(cloneDir_v, x), select = toRead_v)

    ## Get sample (S[0-9]+), sample number ([0-9]+), and treatment
    currSamp_v <- strsplit(x, split = "_")[[1]][2]
    currNum_v <- as.numeric(gsub("S", "", currSamp_v))

    ## Add sample column
    y$Patient <- currSamp_v

    ## Reformat V and J
    y[, (vj_v) := lapply(.SD, function(z) gsub("\\*00|", "", z)), .SDcols = vj_v]
    
    ## Change column order
    if (!is.null(metaFile_v)) {

        ## Get treatment
        currTreat_v <- meta_dt[get(sampleCol_v) == currNum_v, get(treatCol_v)]

        ## Add to data.table
        y$Treatment <- currTreat_v

        ## Re-order columns
        y <- y[,c(1:3,5,4,6), with = F]

        ## Rename columns
        colnames(y) <- c("CDR3b", "TRBV", "TRBJ", "Patient", "Counts", "Treatment")

    } else {
        ## Re-order columns
        y <- y[,c(1:3,5,4), with = F]
        ## Change column names
        colnames(y) <- c("CDR3b", "TRBV", "TRBJ", "Patient", "Counts")
    } # fi

    ## Return
    return(y)

}, simplify = F)

### Combine Together
clones_dt <- do.call("rbind", clones_lsdt)

if (is.null(metaFile_v)) {
    ### Write Output
    outName_v <- paste0(batchName_v, "_gliphClones.txt")
    outFile_v <- file.path(outDir_v, outName_v)

    write.table(clones_dt, outFile_v, row.names = F, quote = F, sep = '\t')
} else {
    ### Get treatments
    treatments_v <- unique(meta_dt[,get(treatCol_v)])
    ### Write one file for each
    for (i in 1:length(treatments_v)){
        currTreat_v <- treatments_v[i]
        currData_dt <- clones_dt[get(treatCol_v) == currTreat_v,]
	currOutName_v <- paste0(batchName_v, "_", currTreat_v, "_gliphClones.txt")
	currOutFile_v <- file.path(outDir_v, currOutName_v)
	write.table(currData_dt[,1:(ncol(currData_dt)-1), with = F], currOutFile_v, row.names = F, quote = F, sep = '\t')
    } # for
} # fi

#for (i in 1:length(cloneFiles_v)){
#
#  ### Get data and names
#  currData_dt <- fread(file.path(cloneDir_v, cloneFiles_v[i]))
#  currName_v <- strsplit(cloneFiles_v[i], split = "_"); currBatch_v <- currName_v[[1]][1]; currSample_v <- currName_v[[1]][2]
#
#  ### Reformat V and J columns
#  columns_v <- c("Best V hit", "Best J hit")
#  currData_dt[, (columns_v) := lapply(.SD, function(x) gsub("\\*00", "", x)), .SDcols = columns_v]
#
#  ### Get V, J, and CDR3
#  currOut_dt <- data.table(currData_dt$`AA. Seq. CDR3`, currData_dt$`Best V hit`, currData_dt$`Best J hit`)
#
#  ### Add Column Name
#  colnames(currOut_dt) <- c("CDR3b", "TRBV", "TRBJ")
#
#  ### Output Names
#  currOutName_v <- paste(currBatch_v, currSample_v, "gliphClones.txt", sep = "_")
#  
#  ### Write Output
#  write.table(currOut_dt, file.path(outDir_v, currOutName_v), row.names = F, quote = F, sep = '\t')
#}
