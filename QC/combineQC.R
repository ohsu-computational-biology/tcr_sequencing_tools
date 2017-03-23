###
### Combine all of the QC files into an excel workbook
###

library(data.table)
library(xlsx)

### Arguments
arguments <- commandArgs(trailingOnly = T)

input.dir <- arguments[1]
output.dir <- arguments[2]

### Get files
input.files <- list.files(input.dir)

### Read and add to workbook
for (i in 1:length(input.files)){
    ## Get file and name
    curr.file <- input.files[i]
    curr.name <- strsplit(curr.file, split = "\\.|_")[[1]]
    print(curr.name)

    ## Skip non-txt files
    if (curr.name[length(curr.name)] != "txt") next

    ## Read input - pear full log is weird, so just don't read it.
    if (!(curr.name[1] == "pear" & curr.name[2] == "full")){
        curr.data <- fread(file.path(input.dir, curr.file))
    }
    

    ## Determine output name
    if (curr.name[1] == "aggregate") {
        sheet.name <- "normalization"
    } else if (curr.name[1] == "count"){
        sheet.name <- curr.name[3]
    } else if (curr.name[1] == "mixcr"){
        sheet.name <- curr.name[2]
    } else if (curr.name[1] == "pear"){
        if (curr.name[2] == "full") {
            next
        } else {
            sheet.name <- "pear"
        }
    } else if (curr.name[1] == "remove"){
        sheet.name <- "remove"
    } else if (curr.name[1] == "uniques"){
        sheet.name <- "analysis"
    } else if (curr.name[2] == "contaminationQC"){
        sheet.name <- "decontam"
    } else if (curr.name[2] == "treatment"){
        next
    } else if (curr.name[2] == "contam"){
        next
    } else if (curr.name[1] == "nonStandard"){
        sheet.name <- paste("nonStandard", curr.name[2], sep = '_')
    } # fi
        

    ## Determine if create workbook or append to workbook
    if (i == 1){
        append.log = F
    } else {
        append.log = T
    } # fi

    ## write sheet
    write.xlsx(curr.data, file = file.path(output.dir, "QC_and_analysis.xlsx"), sheetName = sheet.name, row.names = F, append = append.log)
    
}
