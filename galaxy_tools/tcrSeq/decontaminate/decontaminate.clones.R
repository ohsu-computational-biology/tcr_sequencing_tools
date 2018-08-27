#!/usr/bin/env Rscript

####################################################
### Remove monoclonal contamination from samples ###
####################################################

suppressMessages(library(data.table))
suppressMessages(library(optparse))

# We have three monoclonal sequences that have been used at various points throughout the project. After their introduction,
# we have noticed an overabundance of their presence in later samples. This script searches through clone files exported by
# Mixcr exportClones using the V identity, J identity, and AA sequence of the three clonal contaminants and removes them from
# the file. The output is the exact same as our usual MiXCR export, minus these clones

##################
### CLONE INFO ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################

contam_lsv <- list("p14" = c("V133", "J2-4", "CASSDAGGRNTLYF"),
                   "ot1" = c("V121", "J2-7", "CASSRANYEQYF"),
                   "el4" = c("V15", "J2-3", "CASSTGTETLYF"))

####################
### COMMAND LINE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

optlist <- list(
  make_option(
    c("-c", "--cloneInputs"),
    type = "character",
    help = "MiXCR export clone files. Commonly un-normalized, but can be normalized."
  ),
  make_option(
    c("-n", "--cloneNames"),
    type = "character",
    help = "Actual file names of inputs (rather than dataset_000.dat"
  ),
  make_option(
    c("-s", "--sampleID"),
    type = "integer",
    help = "Index of sample number to extract from sample name, when name is split by '_'."
  ),
  make_option(
    c("-o", "--cloneOutput"),
    type = "character",
    help = "MiXCR clone file of same format as cloneInputs, but with contaminated clones removed."
  ),
  make_option(
    c("-q", "--qcOutput"),
    type = "character",
    help = "Tab-separated QC file detailing which clones were removed in each file and their abundances (both absolute and relative)."
  )
)

### Parse Command Line
p <- OptionParser(usage = "%prog -c cloneInputs -n cloneNames -s sampleID -o cloneOutput -q qcOutput",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

### Get command-line arguments
clone_inputs <- args$cloneInputs
clone_name <- args$cloneNames
sample_id <- args$sampleID
clone_outputs <- args$cloneOutput
qc_output <- args$qcOutput

### Empty QC Matrix
contamination.qc <- matrix(nrow = 1, ncol = 11)

###############
### WRANGLE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###############

### Turn id index into numeric
sample_id <- as.numeric(sample_id)

### Extract sample number from name
sample_number <- strsplit(clone_name, split = "_")[[1]][sample_id]

### Read in data
curr.clone <- fread(clone_inputs)

### Remove duplicate columns
colCounts_dt <- as.data.table(table(colnames(curr.clone)))
dupCols_dt <- colCounts_dt[N > 1,]
removeIndexes_v <- unlist(sapply(dupCols_dt$V1, function(x) {
  y <- grep(x, colnames(curr.clone))
  z <- y[2:length(y)]
  return(z)}, simplify = F))
curr.clone <- curr.clone[,-removeIndexes_v,with=F]

### Grab columns
idCol_v <- grep("[Ii][dD]", grep("[Cc]lone", colnames(curr.clone), value = T), value = T)
countCol_v <- grep("[Nn]orm|nb", grep("[Cc]ount", grep("[Cc]lone", colnames(curr.clone), value = T), value = T), value = T, invert = T)
freqCol_v <- grep("[Ff]raction", grep("[Cc]lone", colnames(curr.clone), value = T), value = T)
hitCols_v <- grep("[Hh]it", grep("[Bb]est", colnames(curr.clone), value = T), value = T)
vHit_v <- grep("V", hitCols_v, value = T); jHit_v <- grep("J", hitCols_v, value = T)
newV_v <- "V segments"; newJ_v <- "J segments"
aaSeqCol_v <- grep("[Aa]{2}", grep("CDR3", colnames(curr.clone), value = T), value = T)
readsCol_v <- grep("[Rr]eads", colnames(curr.clone), value = T)

### Add index for ranking
curr.clone$index <- seq(1, length(curr.clone[[idCol_v]]))

## Get first QC info
orig.unique.count <- length(curr.clone[[idCol_v]])
orig.total.count <- sum(curr.clone[[countCol_v]])

## Update column names and contents for V and J segments
curr.clone[[newV_v]] <- gsub("TRB|\\*00", "", curr.clone[[vHit_v]])
curr.clone[[newV_v]] <- gsub("-", "", curr.clone[[newV_v]])
curr.clone[[newJ_v]] <- gsub("TRB|\\*00", "", curr.clone[[jHit_v]])

###################
### OUTPUT PREP ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################

### Counter for how many reads are devoted to contaminated sequences
totalContamCount_v <- 0

### Empty data.table to contain results
qcOut_dt <- data.table("Sample" = sample_number, "Orig.Uniq.Clones" = orig.unique.count, "Orig.Total.Clones" = orig.total.count, "Remaining.Clones" = numeric(length = 1),
                       "Contam.Clones" = numeric(length = 1))
contamCols_v <- unlist(sapply(names(contam_lsv), function(x) paste(x, c("rank", "count"), sep = '.'), simplify = F))
for (col_v in contamCols_v) qcOut_dt[[col_v]] <- numeric(length = 1)

### Copy of clone data for removal
### TODO - can I do this without making a copy?
forRemoval_dt <- curr.clone

################
### DECONTAM ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################

### For each previously-identified contaminant, extract any matching entries
### Take just the first entry, if more than one is found
### Record its rank and count, remove from data, and update rolling count of total contaminant clones

for (i in 1:length(contam_lsv)){
  ## Get current info
  currContam_v <- names(contam_lsv)[i]; print(sprintf("Currently checking for: %s", currContam_v))
  currV_v <- contam_lsv[[currContam_v]][1]
  currJ_v <- contam_lsv[[currContam_v]][2]
  currCDR3_v <- contam_lsv[[currContam_v]][3]
  
  ## Get appropriate output columns
  currRankCol_v <- grep(currContam_v, grep("rank", colnames(qcOut_dt), value = T), value = T)
  currCountCol_v <- grep(currContam_v, grep("count", colnames(qcOut_dt), value = T), value = T)
  
  ## Extract clone
  currOffendingClone_dt <- curr.clone[(get(newV_v) == currV_v &
                                         get(newJ_v) == currJ_v &
                                         get(aaSeqCol_v) == currCDR3_v),]
  
  ## Subset for first only
  currOffendingClone_dt <- currOffendingClone_dt[1,]
  
  ## Update count
  if (is.na(currOffendingClone_dt[[readsCol_v]][1])) {
    currContamCount_v <- 0
  } else {
    currContamCount_v <- length(unlist(strsplit(as.character(currOffendingClone_dt[[readsCol_v]]), split = ",")))
  } # fi
  
  ## Update total
  totalContamCount_v <- totalContamCount_v + currContamCount_v
  
  ## Get rank
  currRank_v <- currOffendingClone_dt$index[1]
  
  ## Update rank and count
  qcOut_dt[1,eval(currRankCol_v)] <- currRank_v
  qcOut_dt[1,eval(currCountCol_v)] <- currContamCount_v
  
  ## Add contaminant entry back to data.table for later removal
  forRemoval_dt <- rbind(forRemoval_dt, currOffendingClone_dt)
}

### Remove clones
final.clone <- forRemoval_dt[!duplicated(forRemoval_dt, fromLast = T) & seq(nrow(forRemoval_dt)) <= nrow(curr.clone),]
final.clone <- final.clone[!(is.na(final.clone[[idCol_v]]))]

### Final QC update
remainingCount_v <- orig.total.count - totalContamCount_v

qcOut_dt$Remaining.Clones <- remainingCount_v
qcOut_dt$Contam.Clones <- totalContamCount_v

#######################
### POST-PROCESSING ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################

## Re-calculate the clone fractions of the removed clones

new_total_count_v <- sum(final.clone[[countCol_v]])
final.clone[[freqCol_v]] <- final.clone[[countCol_v]] / new_total_count_v

##############
### OUTPUT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############

### Write decontaminated clone
write.table(final.clone, clone_outputs,
            sep = '\t', quote = F, row.names = F)

### Write QC row
write.table(contamination.qc, qc_output,
            sep = '\t', quote = F, row.names = F)
