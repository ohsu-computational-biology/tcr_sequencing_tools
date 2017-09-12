###
### Convert MiXCR to GLIPH
###

### Convert clone files in the MiXCR format into GLIPH-compatible tables

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
