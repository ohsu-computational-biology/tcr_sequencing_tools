###
### CREATE GLIPH INPUT
###

###################
### DESCRIPTION ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################

### For a given batch, first create an aggregate file of all normalized clones and their frequency group designations
### Must also have treatment designations as well.
### For each treatment, create a new file with maximum 10k unique clones in the GLIPH input format
### Determine which 10k clones to include by starting with most frequent to least frequent (by frequency group)
### Within a group, prioritize higher degree of sharing among samples.

####################
### DEPENDENCIES ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

library(wrh.rUtils)
library(optparse)

####################
### COMMAND LINE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

### Make list of options
optlist <- list(
  make_option(
    c("-i", "--inputFile"),
    type = "character",
    help = "'full' clone file from frequency group output"),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    help = "path to directory to write output files"),
  make_option(
    c("-n", "--nClones"),
    type = "numeric",
    default = 10000,
    help = "Maximum number of unique clones allowed in each treatment group"),
  make_option(
    c("-f", "--freqGroups"),
    type = "character",
    default = "Hyperexpanded,Large,Medium,Small,Rare",
    help = "list of frequency groups to use. This cut-off takes precedence over nClones cut-off. Comma-sep, no spaces"
  )
)

### Parse commandline
p <- OptionParser(usage = "%prog -i inputFile -o outDir -n nClones -f freqGroups",
                  option_list = optlist)

args <- parse_args(p)
opt <- args$options

### Assign
inputFile_v <- args$inputFile
outDir_v <- args$outDir
nClones_v <- args$nClones
freqs_v <- args$freqGroups

#############
### SETUP ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############

### Get data
inputData_dt <- fread(inputFile_v)

### Get treatments
treats_v <- unique(inputData_dt$Treatment)

### Frequency groups
freqs_v <- splitComma(freqs_v)

### Output
out_lsdt <- list()

###########
### RUN ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########

for (i in 1:length(treats_v)) {
  
  ## Get treatment and subset data
  currTreat_v <- treats_v[i]
  currData_dt <- inputData_dt[Treatment == currTreat_v,]
  cat(sprintf("Currently on %s (%s of %s)\n", currTreat_v, i, length(treats_v)))
  
  ## Empty data table
  currOut_dt <- data.table()
  stop_v <- F
  
  ## Run for each freq
  for (j in 1:length(freqs_v)) {
    
    ## Initial check
    if (stop_v) next
    
    ## Get freq and subset
    currFreq_v <- freqs_v[j]
    currFreqData_dt <- currData_dt[Div == currFreq_v,]
    
    ## Skip if empty
    if (currFreqData_dt[,.N] == 0) next
    
    ## Create table of occurrences
    currTable_dt <- as.data.table(table(currFreqData_dt[,mget(c("aaSeqCDR3", "V segments", "J segments"))]))[N > 0,][rev(order(N))]
    
    ## Add
    for (k in 1:currTable_dt[,.N]) {
      
      ## Initial check
      if (stop_v) next
      
      ## Get clones
      currAdd_dt <- currFreqData_dt[aaSeqCDR3 == currTable_dt[k,aaSeqCDR3] &
                                      `V segments` == currTable_dt[k,`V segments`] &
                                      `J segments` == currTable_dt[k,`J segments`],
                                    mget(c("aaSeqCDR3", "V segments", "J segments", "Sample", "nb.clone.count"))]
      
      ## Change sample
      # currAdd_dt$Sample <- gsub("^S", "", currAdd_dt$Sample)
      # currAdd_dt$Sample <- sapply(currAdd_dt$Sample, function(x) paste0(paste0(rep(0, (3-nchar(x))), collapse = ""), x))
      # 
      # ## Change count
      # currAdd_dt$nb.clone.count <- sapply(currAdd_dt$nb.clone.count, function(x) paste0(paste0(rep(0, (7-nchar(x))), collapse = ""), x))
      # 
      # ## Make new patientCounts column
      # currAdd_dt$patientCounts <- paste(currAdd_dt$Sample, currAdd_dt$nb.clone.count, sep = "/")
      # 
      # ## Remove other columns
      # rmCol_v <- c("nb.clone.count", "Sample")
      # currAdd_dt[, (rmCol_v) := NULL]
      
      ## Change column names
      colnames(currAdd_dt) <- c("CDR3b", "TRBV", "TRBJ", "Patient", "Counts")
      
      if ((nrow(currOut_dt) + nrow(currAdd_dt)) < nClones_v) {
        currOut_dt <- rbind(currOut_dt, currAdd_dt)
      } else {
        stop_v <- T
        cat(sprintf("Maximum reached. Stopped at %s clone of %s freq group.\n", k, currFreq_v))
        #next
      } # fi
      
    } # for k
    
  } # for j
  
  ## Add to list
  out_lsdt[[currTreat_v]] <- currOut_dt
  
} # for i

### Write
sapply(names(out_lsdt), function(x) write.table(out_lsdt[[x]], file = file.path(outDir_v, paste0(x, "_gliphClones.txt")), 
                                                sep = '\t', quote = F, row.names = F))