###
### Create GLIPH reference library
###

### Using wild-type "naive" TCR samples, create a reference set for use in the GLIPH algorithm

### Output is a FASTA file where the header line has V, J, and CDR3 information and the sequence line is the CDR3 sequence
### Each TCR in a sample will be listed 1 time

### Header line:
### Starts with ">"
### Column 1 : VNAME, JNAME, CDR3 AA SEQ
### Column 2 : VNAME 300 0 
  ### format : TRBV12-1
  ### input column : Best V hit
### Column 3 : blank
### Column 4 : JNAME 30 0 
  ### format : TRBJ2-5
  ### input column : Best J hit
### Column 5 : CDR3 AA SEQ
  ### input column : AA. Seq. CDR3
### Semicolon separated
### Sequence line:
### CDR3 AA SEQ

### Make list of options
library(optparse)
optlist <- list(
  make_option(
    c("-i", "--inputDir"),
    type = "character",
    default = "~/src/gliph/customDB/inputSeqs",
    help = "Directory of MiXCR clone files"
  ),
  make_option(
    c("-o", "--outFile"),
    type = "character",
    default = "~/src/gliph/db/murineMammaryNaive.fa",
    help = "file path and name of output file"
  )
)

### Parse commandline
p <- OptionParser(usage = "%prog -i inputDirectory -o outputFile",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

### Dependencies
suppressMessages(library(data.table))
suppressWarnings(suppressMessages(library(Biostrings, warn.conflicts = F)))

### Commands
cloneDir_v <- args$inputDir
outFile_v <- args$outFile

### Get files
cloneFiles_v <- list.files(cloneDir_v)

### Initialize output data.table
overall_dt <- NULL

for (i in 1:length(cloneFiles_v)){
  
  ### Get data
  currData_dt <- fread(file.path(cloneDir_v, cloneFiles_v[i]))
  
  ### Reformat V and J columns
  columns_v <- c("Best V hit", "Best J hit")
  currData_dt[, (columns_v) := lapply(.SD, function(x) gsub("\\*00", " 300 0", x)), .SDcols = columns_v]
  currData_dt$`Best J hit` <- gsub("300", "30", currData_dt$`Best J hit`)
  
  ### Get V, J, and CDR3
  currTemp_dt <- data.table(currData_dt$`Best V hit`, "", currData_dt$`Best J hit`, currData_dt$`AA. Seq. CDR3`)
  
  ### Combine
  overall_dt <- rbind(overall_dt, currTemp_dt)
}

### Get unique sequences only
setkey(overall_dt, NULL)
overall_dt <- unique(overall_dt)

### Combine for first column
combo_v <- paste(strsplit(overall_dt$V1, split = ' ')[[1]][1], strsplit(overall_dt$V3, split = ' ')[[1]][1], overall_dt$V4, sep = ',')

### Combine everything into header
header_v <- paste(combo_v, overall_dt$V1, overall_dt$V2, overall_dt$V3, overall_dt$V4, sep = ';')

### Named Vector of sequences
seqs_v <- overall_dt$V4; names(seqs_v) <- header_v

### Turn into string set
StringSet <- BStringSet(seqs_v)

### Write out
Biostrings::writeXStringSet(StringSet, filepath = outFile_v)
  
