#!/usr/bin/Rscript

###
### Convert MiXCR to GLIPH
###

### We can take MiXCR files and then group them by their different clonal frequencies. Sometimes we want to then put these clones
### through the GLIPH algorithm. These files are significantly different in format than the original MiXCR files, so they require a 
### separate script in order to convert them into GLIPH-compatible tables. 
### Since these files have already been grouped by treatment, there will be one input directory and then one output file for each file in that directory.

### Format:
  ### Column 1 : CDR3 AA seq
    ### freqGrp - "AA. Seq. CDR3"
    ### GLIPH - "CDR3b"
  ### Column 2 : V sequence
    ### freqGrp - "V segments"; format V131
    ### GLIPH - "TRBV"; format TRBV13-1
  ### Column 3 : J sequence
    ### freqGrp - "J segments"; format J2-5
    ### GLIPH - "TRBJ"; format TRBJ2-5
  ### Column 4 : Patient (Sample, but want consistent names with tool)
    ### freqGrp - "Sample"; format S1
  ### Column 5 : Counts
    ### freqGrp - "Normalized clone count"


### Dependencies
library(data.table)
library(optparse)

### Make list of options
optlist <- list(
  make_option(
    c("-i", "--inputDir"),
    type = "character",
    help = "Directory of freqGroup clone files. One file for each treatment, generally."
  ),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    help = "output directory to write GLIPH-compatible clone files"
  )
)

### Parse commandline
p <- OptionParser(usage = "%prog -i inputDirectory -o outputDirectory",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

### Commands
cloneDir_v <- args$inputDir
outDir_v <- args$outDir

### Get files
cloneFiles_v <- list.files(cloneDir_v)

### Get batch
batchName_v <- strsplit(cloneFiles_v[1], split = "_")[[1]][1]

### Column variables
toRead_v <- c("AA. Seq. CDR3", "V segments", "J segments", "Normalized clone count", "Sample")
vj_v <- c("V segments", "J segments")
weirdV_v <- c("V121", "V122", "V131", "V132", "V133")

### Read in files and fix up
clones_lsdt <- sapply(cloneFiles_v, function(x) {
    ## Get data
    y <- fread(file.path(cloneDir_v, x), select = toRead_v)

    ## Change sample column to Patient
    colnames(y)[colnames(y) == "Sample"] <- "Patient"

    ## Add dash to specific V's
    y[, `V segments` := gsub('^(V1[23]{1})([123]{1})$', '\\1-\\2', `V segments`)]
    
    ## Reformat V and J
    y[, (vj_v) := lapply(.SD, function(z) paste0("TRB", z)), .SDcols = vj_v]

    ## Re-order columns
    y <- y[,c(1:3,5,4), with = F]
    ## Change column names
    colnames(y) <- c("CDR3b", "TRBV", "TRBJ", "Patient", "Counts")

    ## Return
    return(y)

}, simplify = F)

for (i in 1:length(clones_lsdt)){
    ## Get output name
    currTreat_v <- strsplit(names(clones_lsdt)[i], split = "_")[[1]][2]
    currOut_v <- paste0(batchName_v, "_", currTreat_v, "_gliphClones.txt")
    ## Write table
    write.table(clones_lsdt[[i]], file.path(outDir_v, currOut_v), row.names = F, quote = F, sep = '\t')
} # for 
