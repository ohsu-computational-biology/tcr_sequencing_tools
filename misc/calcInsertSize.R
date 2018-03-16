### Calculate insert size of PEAR-assembled files
### Input is a directory of files that have the read lengths of all the sequences in a PEAR-assembled file
### Each file should be made via `awk '{if (NR % 4 == 2) print length($0)}' <file>`

library(data.table)

arguments <- commandArgs(trailingOnly = T)
inputDir_v <- arguments[1]
outDir_v <- arguments[2]

### Get files
inputFiles_v <- list.files(inputDir_v, pattern = "*.txt")
inputFiles_v <- inputFiles_v[order(as.numeric(gsub(".*_S|.txt", "", inputFiles_v)))]

### Create summary matrix
summary_mat <- matrix(nrow = length(inputFiles_v), ncol = 2)

### Count insert size for each file
for (i in 1:length(inputFiles_v)){
	## Get file and data
	currFile_v <- inputFiles_v[i]
	currData_dt <- fread(file.path(inputDir_v, currFile_v))
	## Calc average and SD
	currAvg_v <- round(mean(currData_dt$V1), digits = 1)
	currSD_v <- round(sd(currData_dt$V1), digits = 1)
	## Update matrix
	summary_mat[i,] <- c(currAvg_v, currSD_v)
} # for i

### Change to data.table and add names
summary_dt <- as.data.table(summary_mat)
summary_dt$File <- inputFiles_v
summary_dt <- summary_dt[,c(3,1,2),with=F]
colnames(summary_dt) <- c("File", "AvgInsert", "SDInsert")

### Write
write.table(summary_dt, file.path(outDir_v, "insertSummary.txt"), row.names = F, quote = F, sep = '\t')
