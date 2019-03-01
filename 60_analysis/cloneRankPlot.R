#!/ usr/bin/Rscript

###
### PLOT CLONE RANKS 
###

###################
### DESCRIPTION ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################

### Sort clones in a sample in decreasing order and plot them with rank along x-axis and count along y-axis.

### Most basic version is 1 plot per sample.

### Second version is 1 plot per mouse, which means multiple lines on the plot for each sample from that mouse.

####################
### DEPENDENCIES ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

library(wrh.rUtils)
loadLib(c("data.table", "optparse", "RColorBrewer", "VennDiagram", "MASS", "scales", "xlsx", "ggplot2"))
#options(scipen=999) # default is 0
options(scipen=0)

####################
### COMMAND LINE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

### Make list of options
optlist <- list(
  make_option(
    c("-i", "--inputFile"),
    type = "character",
    help = "full_clones.txt file created by groupClones.R"
  ),
  make_option(
    c("-m", "--meta"),
    type = "character",
    help = "Metadata file containing treatment designations at a minimum. Require: Sample | Mouse/Animal | Division 2 | Division 1 (optional)"
  ),
  make_option(
    c("-p", "--pairCol"),
    type = "character",
    help = "Metadata column name used to check correct pairing. Could be sample number if separate batches, 
            or mouse/animal ID if same batch. Comma-sep, no spaces if more than one."
  ),
  make_option(
    c("-P", "--whichPairs"),
    type = "character",
    help = "Comma-separated list of values from the 'pairs' column in metadata. This is used to select which samples to compare.
    Example: if 'pairs' is the Mouse column, then '123,124,125' would select the three mice with those IDs to run through the analysis.
    NULL will use all."
  ),
  make_option(
    c("-c", "--compareCol"),
    type = "character",
    default = "tissue",
    help = "Metadata column used to compare overlaps. Each mouse needs to have multiple different values for this. 
          'tissue' (default) will compare by tissue. 'time' is used for longitudinal. This is Division 2."
  ),
  make_option(
    c("-C", "--whichCompare"),
    type = "character",
    help = "Comma-separated list of values from the 'compare' column in metadata. This is used to select which groups to compare.
    Example: if 'compare' column is tissue, then 'tumor,blood' will only select those tissues to compare, even if data exists for 'lung' also.
    Example: if 'compare' column is time, then '1m,endpoint' will only select those time points to compare, even if data exists for '4m' also.
    NULL will use all."
  ),
  make_option(
    c("-s", "--sepCol"),
    type = "character",
    help = "Name of column to use to differentiate between groups of samples. Common to use 'treatment' here. this is Division 1."
  ),
  make_option(
    c("-M", "--merge"),
    type = "character",
    help = "'no' - one plot per mouse; 'merge' - merge all mice together and create new ranks; 'mean' - merge mice together by taking mean of each rank."
  ),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    help = "file path to directory for writing output files"
  ),
  make_option(
    c("-d", "--debug"),
    type = "logical",
    default = FALSE,
    help = "Logical. TRUE - print session info and extra output to help with debugging. 
            Also do not write output (tables and images). FALSE - normal output and write output files (tables and images)."
  ),
  make_option(
    c("-l", "--log"),
    type = "logical",
    default = FALSE,
    help = "Logical. TRUE - output session info. FALSE - do not output session info. 
          If debug is TRUE and log is FALSE, then session info will be printed to STDOUT. If neither are set, no output."
  ),
  make_option(
    c("-w", "--write"),
    type = "logical",
    default = TRUE,
    help = "logical indicating whether or not to save plots."
  )
)



### Parse commandline
p <- OptionParser(usage = "%prog -i inputDir -m meta -p pairCol -P whichPairs -c compareCol -C whichCompare -s sepCol 
                  -M merge -o outDir -d debug -l log -w write",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

### Assign arguments
inputFile_v <- args$inputFile
metaFile_v <- args$meta
pairCol_v <- args$pairCol
whichPair_v <- args$whichPairs
compareCol_v <- args$compareCol
whichCompare_v <- args$whichCompare
sepCol_v <- args$sepCol
merge_v <- args$merge
outDir_v <- args$outDir
debug_v <- args$debug
log_v <- args$log
toWrite_v <- args$write

### Handle which arguments
if (!is.null(whichPair_v)) whichPair_v <- splitComma(whichPair_v)
if (!is.null(whichCompare_v)) whichCompare_v <- splitComma(whichCompare_v)

### Check arguments
print("Current input file:")
print(inputFile_v)

#############
### SETUP ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############

### Print commands
if (log_v){
  returnSessionInfo(args_lsv = args, out_dir_v = outDir_v)
} else {
  if (debug_v){
    returnSessionInfo(args_lsv = args)
  } # fi
} # fi

### Get data
meta_dt <- fread(metaFile_v)
inputData_dt <- fread(inputFile_v)

### Get column names
sampleCol_v <- grep("[Ss]ample", colnames(meta_dt), value = T)

## Get pair column name
## This is the column that identifies which mice are pairs
pairCol_v <- unlist(sapply(pairCol_v, function(x){
  y <- grep(x, colnames(meta_dt), value = T)
  if (length(y) == 0) {
    y <- grep(simpleCap(pairCol_v), colnames(meta_dt), value = T)
  }
  return(y)
}, simplify = F))

### Get clone identity columns
seqCol_v <- grep("nSeqCDR3|clonalSequence", colnames(inputData_dt), value = T)
vCol_v <- grep("V Segments|V segments", colnames(inputData_dt), value = T)
jCol_v <- grep("J Segments|J segments", colnames(inputData_dt), value = T)

### Get clone count column
countCol_v <- "nb.clone.count"

### Add compare column
inputData_dt$compare <- paste(inputData_dt[[vCol_v]], inputData_dt[[seqCol_v]], inputData_dt[[jCol_v]], sep = "_")

### Make sure compareCol_v is in inputdata
if (!compareCol_v %in% colnames(inputData_dt)) {
  inputData_dt <- merge(inputData_dt, meta_dt[,mget(c(sampleCol_v, compareCol_v))], by = sampleCol_v, sort = F, all = F)
}

####################
### SUBSET INPUT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

### Subset for whichPair and whichCompare
if (!is.null(whichPair_v)) {
  subMeta_dt <- meta_dt[get(pairCol_v) %in% whichPair_v,]
} else {
  subMeta_dt <- meta_dt
}

if (!is.null(whichCompare_v)) {
  subMeta_dt <- subMeta_dt[get(compareCol_v) %in% whichCompare_v,]
} else {
  subMeta_dt <- subMeta_dt
}

####################
### PREPARE DATA ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

sepData_lsdt <- sepMeta_lsdt <- list()

if (is.null(sepCol_v)) {
  
  sepData_lsdt[["all"]] <- inputData_dt
  subMeta_dt$master <- paste(subMeta_dt[[pairCol_v]], subMeta_dt[[sampleCol_v]], subMeta_dt[[compareCol_v]], sep = "_")
  sepMeta_lsdt[["all"]] <- subMeta_dt
  
} else {
  
  ## Get unique values
  sepVals_v <- unique(meta_dt[[sepCol_v]])
  
  ## Subset
  for (val_v in sepVals_v) {
    
    ## Get samples
    currSamples_v <- subMeta_dt[get(sepCol_v) == val_v, get(sampleCol_v)]
    currSamples_v <- paste0("S", gsub("^S", "", currSamples_v))
    
    ## Subset metadata and add column
    currMeta_dt <- subMeta_dt[get(sampleCol_v) %in% currSamples_v,]
    currMeta_dt$master <- paste(currMeta_dt[[pairCol_v]], currMeta_dt[[sampleCol_v]], 
                                currMeta_dt[[compareCol_v]], currMeta_dt[[sepCol_v]], sep = "_")
    
    ## Subset data
    sepData_lsdt[[val_v]] <- inputData_dt[Sample %in% currSamples_v,]
    sepMeta_lsdt[[val_v]] <- currMeta_dt
    
  } # for val_v
} # fi

############################
### PERFORM CALCULATIONS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############################

for (i in 1:length(sepData_lsdt)) {
  
  ## Get name, data, and possible comparisons
  currName_v <- names(sepData_lsdt)[i]
  currData_dt <- sepData_lsdt[[currName_v]]
  currMeta_dt <- sepMeta_lsdt[[currName_v]]
  
  ## Get individual animals
  currMice_v <- unique(currMeta_dt[[pairCol_v]])
  
  ## Make sub-directory, if needed
  if (currName_v != "all" & !(merge_v %in% c("merge", "mean"))) {
    currDir_v <- mkdir(outDir_v, currName_v)
  } else {
    currDir_v <- outDir_v
  }
  
  ## Combine, if needed
  if (merge_v == "merge") {
    currMeta_dt$merge <- rep("merge", nrow(currMeta_dt))
    currMice_v <- "merge"
  } else if (merge_v == "mean") {
    currMeta_dt$merge <- rep("mean", nrow(currMeta_dt))
    currMice_v <- "mean"
  }
  
  ## Run for each animal
  for (j in 1:length(currMice_v)) {
    
    ###
    ### GET SPECIFIC MOUSE INFO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###
    
    ## Get mouse
    currMouse_v <- currMice_v[j]
    
    ## Subset metadata and get samples
    if (merge_v %in% c("merge", "mean")) {
      currMouseMeta_dt <- currMeta_dt
    } else {
      currMouseMeta_dt <- currMeta_dt[get(pairCol_v) == currMouse_v,]
    }
    
    ### Skip if no pairs
    if (nrow(currMouseMeta_dt) <= 1) {
      cat(sprintf("Skipping mouse %s because no pairs\n", currMouseMeta_dt[[pairCol_v]]))
      next
    }
    
    ###
    ### TURN INTO LIST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###
    
    currMouseData_lsdt <- list()
    
    for (k in 1:currMouseMeta_dt[,.N]) {
      
      ## Get specific mouse infomration
      currMouseSample_v <- currMouseMeta_dt[k,get(sampleCol_v)]
      currMouseCompare_v <- currMouseMeta_dt[k,get(compareCol_v)]
      
      ## Make name including mouse, sample, and division 2
      currMouseName_v <- paste(currMouse_v, currMouseSample_v, currMouseCompare_v, sep = "_")
      
      ## If multiple sepData (i.e. using division 1) add that
      if (currName_v != "all") currMouseName_v <- paste(currMouseName_v, currName_v, sep = "_")
      
      ## Check that it matches 'master'
      if (!merge_v %in% c("merge", "mean")) {
        if (currMouseName_v != currMouseMeta_dt[k,master]) stop("Mouse Name does not match 'master' name. Please check.")
      }
      
      ## Subset for frequency groups 
      currSubData_dt <- currData_dt[Sample == currMouseSample_v &
                                      get(compareCol_v) == currMouseCompare_v,]
      
      ## Sort
      currSubData_dt <- currSubData_dt[order(get(countCol_v), decreasing = T)]
      
      ## Add rank
      currSubData_dt$Rank <- 1:nrow(currSubData_dt)
      
      ## Add log count
      currSubData_dt$log10Count <- log10(currSubData_dt[[countCol_v]])
      
      ## Add master identifier
      currSubData_dt$master <- rep(currMouseName_v, nrow(currSubData_dt))
      
      ## Add to list
      currMouseData_lsdt[[currMouseName_v]] <- currSubData_dt
      
    } # for k
    
    ##
    ## FINAL WRANGLE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##
    
    ### Combine
    currMouseData_dt <- do.call(rbind, currMouseData_lsdt)
    
    if (merge_v == "merge") {
      
      newRank_lsdt <- lapply(unique(currMouseData_dt[[compareCol_v]]), function(x) {
        temp_dt <- currMouseData_dt[get(compareCol_v) == x,]
        temp_dt <- temp_dt[order(get(countCol_v), decreasing = T)]
        temp_dt$Rank <- 1:nrow(temp_dt)
        return(temp_dt)
      })
      
      currMouseData_dt <- do.call(rbind, newRank_lsdt)
      
    } else if (merge_v == "mean") {
      
      newRank_lsdt <- lapply(unique(currMouseData_dt[[compareCol_v]]), function(x) {
        temp_dt <- currMouseData_dt[get(compareCol_v) == x,mget(c("Rank", "log10Count", compareCol_v))]
        mean_dt <- temp_dt[, .(Mean = mean(log10Count)), by = Rank]
        mean_dt[[compareCol_v]] <- rep(x, nrow(mean_dt))
        mean_dt <- mean_dt[order(Mean, decreasing = T)]
        mean_dt$Rank <- 1:nrow(mean_dt)
        return(mean_dt)
      })
      
      currMouseData_dt <- do.call(rbind, newRank_lsdt)
      colnames(currMouseData_dt) <- c("Rank", "log10Count", compareCol_v)
      
    } # fi
    
    ##
    ## OUTPUT GRAPH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##
    
    ## Make title and file
    if (merge_v == "merge") {
      currTitle_v <- "Merged"
      currFile_v <- "merged_compareRank.pdf"
    } else if (merge_v == "mean") {
      currTitle_v <- "Mean"
      currFile_v <- "mean_compareRank.pdf"
    } else {
      currTitle_v <- paste0("Mouse - ", currMouse_v)
      currFile_v <- paste0("mouse_", currMouse_v, "_compareRank.pdf")
    }
    
    if (currName_v != "all") {
      currTitle_v <- paste0(currTitle_v, " (", currName_v, ")")
      currFile_v <- paste0(currName_v, "_", currFile_v)
    }
    
    ## Calculate p-value
    currT <- t.test(currMouseData_dt[get(compareCol_v) == whichCompare_v[1], log10Count],
                    currMouseData_dt[get(compareCol_v) == whichCompare_v[2], log10Count])
    currP_v <- formatC(currT$p.value, format = "e", digits = 4)
    
    ## Get annotation parameters
    x_v <- max(currMouseData_dt$Rank) * .8
    y_v <- max(currMouseData_dt$log10Count) * .8
    currAnnotate_v <- paste0("p = ", currP_v)
    
    ### Make plot
    currPlot_gg <- ggplot(currMouseData_dt, aes_string(x = "Rank", y = "log10Count", color = compareCol_v)) + 
      geom_point() +
      #ggtitle(paste0("Mouse: ", currMouse_v)) +
      ggtitle(currTitle_v) +
      big_label() + labs(color = "Timepoint") + 
      annotate(geom = "text", x = x_v, y = y_v, label = currAnnotate_v)
    
    ### Save plot
    if (toWrite_v) {
      ggsave(file.path(currDir_v, currFile_v), currPlot_gg, width = 7, height = 7)
    } else {
      print(currPlot_gg)
    } # fi

  } # for j

} # for i

