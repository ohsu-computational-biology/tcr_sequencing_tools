#!/usr/bin/env Rscript

###
### Combine all of the QC files into an excel workbook
###

suppressMessages(library(data.table))
suppressMessages(library(xlsx))
suppressMessages(library(optparse))

####################
### COMMAND LINE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

optlist <- list(
  make_option(
    c("-c", "--countSpikes"),
    type = "character",
    help = "QC output of count spikes and remove spikes steps."
  ),
  make_option(
    c("-l", "--align"),
    type = "character",
    help = "QC output of MiXCR alignment step."
  ),
  make_option(
    c("-s", "--assemble"),
    type = "character",
    help = "QC output of MiXCR assemble step."
  ),
  make_option(
    c("-d", "--decontam"),
    type = "character",
    help = "Aggregated QC files of decontamination step."
  ),
  make_option(
    c("-n", "--normalize"),
    type = "character",
    help = "QC output of normalization step."
  ),
  make_option(
    c("-a", "--analysis"),
    type = "character",
    help = "Output of diversity analysis metrics from analysis step."
  )
)

### Parse Command Line
p <- OptionParser(usage = "%prog -c countSpikes -l align -s assemble -d decontam -n normalize -a analysis",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

### Get command-line arguments
count_spikes_qc <- args$countSpikes
align_qc <- args$align
assemble_qc <- args$assemble
decontam_qc <- args$decontam
norm_qc <- args$normalize
analysis <- args$analysis

input_files <- c("Spikes" = count_spikes_qc, "Align" = align_qc, "Assemble" = assemble_qc,
                 "Decontam" = decontam_qc, "Normalization" = norm_qc, "Analysis" = analysis)

### Read files and add to excel workbook

for (i in 1:length(input_files)){
    ## Get data
    curr_dt <- fread(input_files[i])
    curr_name <- names(input_files)[i]

    ## Determine if create or append 
    if (i == 1){
        append.log = F
    } else {
        append.log = T
    } # fi

    ## Write Sheet
    write.xlsx(curr_dt,
               file = "./temp.xlsx",
               sheetName = curr_name,
               row.names = FALSE,
               append = append.log)
} # for
