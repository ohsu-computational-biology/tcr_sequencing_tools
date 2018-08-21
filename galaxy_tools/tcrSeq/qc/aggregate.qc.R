#!/usr/bin/Rscript

###
### Combine all of the QC files into an excel workbook
###

suppressMessages(library(data.table))
suppressMessages(library(xlsx))

### Arguments
arguments <- commandArgs(trailingOnly=T)

count_spikes_qc <- arguments[1]
align_qc <- arguments[2]
assemble_qc <- arguments[3]
decontam_qc <- arguments[4]
norm_qc <- arguments[5]
analysis <- arguments[6]
#out_file <- arguments[7]
#out_file <- gsub("dat", "xlsx", out_file)

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
