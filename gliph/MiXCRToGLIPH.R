#!/usr/bin/Rscript

###
### Convert MiXCR to GLIPH
###

### Convert clone files in the MiXCR format into GLIPH-compatible tables
### All samples from a batch will be combined into one file

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
toRead_v <- c("AA. Seq. CDR3", "Best V hit", "Best J hit", "Normalized clone count")
vj_v <- c("Best V hit", "Best J hit")

### Read in files and fix up
clones_lsdt <- sapply(cloneFiles_v, function(x) {
    ## Get data
    y <- fread(file.path(cloneDir_v, x), select = toRead_v)
    ## Add sample column
    y$Patient <- strsplit(x, split = "_")[[1]][2]
    ## Reformat V and J
    y[, (vj_v) := lapply(.SD, function(z) gsub("\\*00|", "", z)), .SDcols = vj_v]
    ## Change column order
    y <- y[,c(1:3,5,4), with = F]
    ## Change column names
    colnames(y) <- c("CDR3b", "TRBV", "TRBJ", "Patient", "Counts")
    ## Return
    return(y)
}, simplify = F)

### Combine Together
clones_dt <- do.call("rbind", clones_lsdt)

### Write Output
outName_v <- paste0(batchName_v, "_gliphClones.txt")
outFile_v <- file.path(outDir_v, outName_v)

write.table(clones_dt, outFile_v, row.names = F, quote = F, sep = '\t')



for (i in 1:length(cloneFiles_v)){

  ### Get data and names
  currData_dt <- fread(file.path(cloneDir_v, cloneFiles_v[i]))
  currName_v <- strsplit(cloneFiles_v[i], split = "_"); currBatch_v <- currName_v[[1]][1]; currSample_v <- currName_v[[1]][2]

  ### Reformat V and J columns
  columns_v <- c("Best V hit", "Best J hit")
  currData_dt[, (columns_v) := lapply(.SD, function(x) gsub("\\*00", "", x)), .SDcols = columns_v]

  ### Get V, J, and CDR3
  currOut_dt <- data.table(currData_dt$`AA. Seq. CDR3`, currData_dt$`Best V hit`, currData_dt$`Best J hit`)

  ### Add Column Name
  colnames(currOut_dt) <- c("CDR3b", "TRBV", "TRBJ")

  ### Output Names
  currOutName_v <- paste(currBatch_v, currSample_v, "gliphClones.txt", sep = "_")
  
  ### Write Output
  write.table(currOut_dt, file.path(outDir_v, currOutName_v), row.names = F, quote = F, sep = '\t')
}
