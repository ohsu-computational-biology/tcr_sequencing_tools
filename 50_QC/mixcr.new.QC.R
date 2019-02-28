#	RScript that aggregates reports from MiXCR's alignment tool
#   
#   At this point the script simply accumulates results, but it'd be easy to add
#       some visualization, analysis, etc. once the data is aggregated

##	Get command-line arguments
library(data.table)

arguments <- commandArgs(trailingOnly=TRUE);

reportDir_v <- arguments[1]
outDir_v <- arguments[3]

## Get files
reportFiles_v <- list.files(reportDir_v)

## Sort files
reportFiles_v <- reportFiles_v[order(as.numeric(gsub("^S|_report.*", "", reportFiles_v)))]

## Get sample numbers
reportNumbers_v <- gsub("_.*", "", reportFiles_v)

## Make data.frame
output_dt <- data.table(Sample = reportNumbers_v)
output_dt$Input.Reads <- numeric()
output_dt$Aligned.Reads <- numeric()
output_dt$pct.Aligned <- numeric()
output_dt$Total.Clones <- numeric()
output_dt$Reads.In.Clones <- numeric()
output_dt$pct.Reads.In.Clones <- numeric()

## Iterate over files and extract info
for (i in 1:length(reportFiles_v)){
   ## Get files
   currReport_v <- reportFiles_v[i]

   ## Get data
   currReport_data <- readLines(file.path(reportDir_v, currReport_v))
  
   ## Get Input reads
   currInputReads_v <- grep("Total sequencing reads", currReport_data, value = T)
   currInputReads_v <- as.numeric(gsub("^.*: ", "", currInputReads_v))

   ## Get Aligned Reads
   currAlignedReads_v <- grep("Successfully aligned reads", currReport_data, value = T)
   currAlignedPct_v <- as.numeric(gsub("^.*\\(|%\\)$", "", currAlignedReads_v))
   currAlignedReads_v <- as.numeric(gsub("^.*: | \\(.*$", "", currAlignedReads_v))

   ## Get Clones
   currClones_v <- grep("Final clonotype count", currReport_data, value = T)
   currClones_v <- as.numeric(gsub(".*: ", "", currClones_v))

   ## Get reads used in clones
   currReadsInClones_v <- grep("Reads used in clonotypes, percent of total", currReport_data, value = T)
   currReadsInClonesPct_v <- as.numeric(gsub("^.*\\(|%\\)$", "", currReadsInClones_v))
   currReadsInClones_v <- as.numeric(gsub(".*: | \\(.*$", "", currReadsInClones_v))

   ## Update data.table
   output_dt[i, "Input.Reads" := currInputReads_v]
   output_dt[i, "Aligned.Reads" := currAlignedReads_v]
   output_dt[i, "pct.Aligned" := currAlignedPct_v]
   output_dt[i, "Total.Clones" := currClones_v]
   output_dt[i, "Reads.In.Clones" := currReadsInClones_v]
   output_dt[i, "pct.Reads.In.Clones" := currReadsInClonesPct_v]
}

### Write output
outName_v <- file.path(outDir_v, "mixcr.QC.summary.txt")
write.table(output_dt, outName_v, quote = F, sep = '\t', row.names = F)
