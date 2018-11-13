#	RScript that aggregates reports from MiXCR's alignment tool
#   
#   At this point the script simply accumulates results, but it'd be easy to add
#       some visualization, analysis, etc. once the data is aggregated

##	Get command-line arguments
library(data.table)

arguments <- commandArgs(trailingOnly=TRUE);

alignDir_v <- arguments[1]
assembleDir_v <- arguments[2]
outDir_v <- arguments[3]

## Get files
alignFiles_v <- list.files(alignDir_v)
assembleFiles_v <- list.files(assembleDir_v)

## Sort files
alignFiles_v <- alignFiles_v[order(as.numeric(gsub("^S|_align.*", "", alignFiles_v)))]
assembleFiles_v <- assembleFiles_v[order(as.numeric(gsub("^S|_assemble.*", "", assembleFiles_v)))]

## Get sample numbers
alignNumbers_v <- gsub("_.*", "", alignFiles_v)
assembleNumbers_v <- gsub("_.*", "", assembleFiles_v)

## Make sure they match
mismatch_v <- which(alignNumbers_v != assembleNumbers_v)
if (length(mismatch_v) > 0) stop("Mismtach report files")

## Make data.frame
output_dt <- data.table(Sample = alignNumbers_v)
output_dt$Input.Reads <- numeric()
output_dt$Aligned.Reads <- numeric()
output_dt$Total.Clones <- numeric()
output_dt$Reads.In.Clones <- numeric()

## Iterate over files and extract info
for (i in 1:length(alignFiles_v)){
   ## Get files
   currAlign_v <- alignFiles_v[i]
   currAssemble_v <- assembleFiles_v[i]

   ## Get data
   currAlign_data <- readLines(file.path(alignDir_v, currAlign_v))
   currAssemble_data <- readLines(file.path(assembleDir_v, currAssemble_v))

   ## Get Input reads
   currInputReads_v <- grep("Total sequencing reads", currAlign_data, value = T)
   currInputReads_v <- as.numeric(gsub("^.*: ", "", currInputReads_v))

   ## Get Aligned Reads
   currAlignedReads_v <- grep("Successfully aligned reads", currAlign_data, value = T)
   currAlignedPct_v <- as.numeric(gsub("^.*\\(|%\\)$", "", currAlignedReads_v))
   currAlignedReads_v <- as.numeric(gsub("^.*: | \\(.*$", "", currAlignedReads_v))

   ## Get Clones
   currClones_v <- grep("Final clonotype count", currAssemble_data, value = T)
   currClones_v <- as.numeric(gsub(".*: ", "", currClones_v))

   ## Get reads used in clones
   currReadsInClones_v <- grep("Reads used in clonotypes, percent of total", currAssemble_data, value = T)
   currReadsInClones_v <- as.numeric(gsub(".*: | \\(.*$", "", currReadsInClones_v))

   ## Update data.table
   output_dt[i, "Input.Reads" := currInputReads_v]
   output_dt[i, "Aligned.Reads" := currAlignedReads_v]
   output_dt[i, "Total.Clones" := currClones_v]
   output_dt[i, "Reads.In.Clones" := currReadsInClones_v]
}

### Write output
outName_v <- file.path(outDir_v, "mixcr.QC.summary.txt")
write.table(output_dt, outName_v, quote = F, sep = '\t', row.names = F)
