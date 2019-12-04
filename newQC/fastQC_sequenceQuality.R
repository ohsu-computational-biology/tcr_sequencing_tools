#!/usr/bin/Rscript

###################
### DESCRIPTION ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################

### Summarize the results of the multiqc output.
### Most data is summarized in the multiqc_data/multiqc_fastqc.txt file
### Currently only take a look at "per_base_sequence_quality"
### Want to see how many samples are flagged and what the distribution of R1/R2 files is.

####################
### DEPENDENCIES ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(optparse))

####################
### COMMAND LINE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

### Make command line
optlist <- list(
  make_option(
    c("-i", "--inputFile"),
    type = "character",
    help = "Path to directory of pear log files. Will only read files that satisfy pear_full_log_[0-9]+.txt"
  ),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    help = "Path to directory to write output files. If this argument is not used, will write to inputDir."
  )
)

### Parse command line
p <- OptionParser(usage = "%prog -i inputFile -o outDir",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

### Get arguments
inputFile_v <- args$inputFile
outDir_v <- args$outDir

#############
### SETUP ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############

### Read data
inputData_dt <- fread(inputFile_v)

### Subset for column
inputData_dt <- inputData_dt[,mget(c("Sample", "per_base_sequence_quality"))]

### Add sample number column
sampleNum_v <- sapply(inputData_dt$Sample, function(x) strsplit(x, split = "_")[[1]][2])
inputData_dt$sampleNum <- sampleNum_v

### Add R1/R2 column
pairColumn_v <- sapply(inputData_dt$Sample, function(x) grep("R[1-2]", strsplit(x, split = "_")[[1]], value = T))
inputData_dt$Pair <- pairColumn_v

############################
### IDENTIFY BAD SAMPLES ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############################

lvls_v <- c("pass", "fail", "warn")
split_lsdt <- lapply(lvls_v, function(x) {
  y <- inputData_dt[per_base_sequence_quality == x,]
  y$Sample <- sapply(y$Sample, function(x) paste(strsplit(x, split = "_")[[1]][1:2], collapse = "_"))
  y1 <- y[Pair == "R1",]; y2 <- y[Pair == "R2",]
  z <- merge(y1, y2,
             by = c("Sample", "per_base_sequence_quality", "sampleNum"), sort = F, suffixes = c("_1", "_2"), all = T)
  z$Pair <- paste(z$Pair_1, z$Pair_2, sep = ";")
  rm_v <- c("Pair_1", "Pair_2", "per_base_sequence_quality")
  z[,(rm_v) := NULL]
  z$Pair <- gsub("NA;|;NA", "", z$Pair)
  z <- z[order(as.numeric(gsub("^S", "", sampleNum)))]
  return(z)
  })
names(split_lsdt) <- lvls_v

############################
### CALC RESULTS SUMMARY ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############################

### TABLE OF ALL SAMPLES
fullResults_dt <- as.data.table(table(inputData_dt$per_base_sequence_quality))

### TABLE BY R1/R2
pairResults_dt <- as.data.table(table(inputData_dt[,mget(c("Pair", "per_base_sequence_quality"))]))

### Combine results
outResults_dt <- merge(fullResults_dt, pairResults_dt[Pair == "R1", mget(c("per_base_sequence_quality", "N"))], by.x = "V1", by.y = "per_base_sequence_quality")
outResults_dt <- merge(outResults_dt, pairResults_dt[Pair == "R2", mget(c("per_base_sequence_quality", "N"))], by.x = "V1", by.y = "per_base_sequence_quality")
colnames(outResults_dt) <- c("Class", "All", "R1", "R2")

### Sort
outResults_dt <- outResults_dt[match(c("fail", "warn", "pass"), outResults_dt$Class),]

### Combine
split_lsdt[["summary"]] <- outResults_dt
split_lsdt <- split_lsdt[c("summary", lvls_v)]

### Write out
writexl::write_xlsx(split_lsdt,
                    path = file.path(outDir_v, "sequenceQualitySummary.txt"))
