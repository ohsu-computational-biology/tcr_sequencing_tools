#!/ usr/bin/Rscript

###
### TRACK CLONES 
###

###################
### DESCRIPTION ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################

### Compare the frequency changes of shared clones between mice.

### Use nucleotide sequence to determine similarity.

### DIVISION 1
### Mice are split into different groups based on some descriptor, such as treatment. It is not possible
### For a single mouse to have membership in multiple of these descriptors. 
### A mice must be all one treatment or all the other (for obvious reasons).
### This descriptor is usually used for summarization - how does the average Jaccard for Treat1 compare to Treat2?

### DIVISION 2
### Each mouse will, however, have multiple instances of some other descriptor, such as tissue or time.
### These are what are used to determine the overlap. For example, what is the overlap between TUMOR, LYMPH NODE, and BLOOD
### of a particular mouse. These are made at the individual mouse level, but can also be aggregated by the division described above.

### For each mouse, make a plot where x-axis is Division 2, y-axis is frequency, each point is a clone, with lines between shared clones.

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
    c("-p", "--pairs"),
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
    c("-c", "--compare"),
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
    c("-s", "--separate"),
    type = "character",
    help = "Name of column to use to differentiate between groups of samples. Common to use 'treatment' here. this is Division 1."
  ),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    help = "file path to directory for writing output files"
  ),
  make_option(
    c("-n", "--wkbkName"),
    type = "character",
    help = "Name of output excel workbook."
  ),
  make_option(
    c("-S", "--sheetName"),
    type = "character",
    help = "Name of workbook sheet."
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
    type = "character",
    default = "all",
    help = "Any of 'all', 'summary', or 'plot'"
  )
)



### Parse commandline
p <- OptionParser(usage = "%prog -i inputDir -m meta -p pairs -P whichPairs -c compare -C whichCompare 
                  -s separate -o outDir -n wkbkName -S sheetName -d debug -l log -w write",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

### Assign arguments
inputFile_v <- args$inputFile
metaFile_v <- args$meta
pair_v <- args$pairs
whichPair_v <- args$whichPairs
compare_v <-args$compare
whichCompare_v <- args$whichCompare
sepCol_v <- args$separate
outDir_v <- args$outDir
wkbkName_v <- args$wkbkName
sheetName_v <- args$sheetName
debug_v <- args$debug
log_v <- args$log
toWrite_v <- args$write

### Handle sheet argument
sheetName_v <- ifelse(is.null(sheetName_v), "Sheet1", sheetName_v)

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
animalCol_v <- grep("[Mm]ouse|[Aa]nimal", colnames(meta_dt), value = T)

## Get pair column name
## This is the column that identifies which mice are pairs
pairCol_v <- unlist(sapply(pair_v, function(x){
  y <- grep(x, colnames(meta_dt), value = T)
  if (length(y) == 0) {
    y <- grep(simpleCap(pair_v), colnames(meta_dt), value = T)
  }
  return(y)
}, simplify = F))

### Get clone identity columns
seqCol_v <- grep("nSeqCDR3|clonalSequence", colnames(inputData_dt), value = T)
vCol_v <- grep("V Segments|V segments", colnames(inputData_dt), value = T)
jCol_v <- grep("J Segments|J segments", colnames(inputData_dt), value = T)

### Add compare column
inputData_dt$compare <- paste(inputData_dt[[vCol_v]], inputData_dt[[seqCol_v]], inputData_dt[[jCol_v]], sep = "_")

### Make output files
allFile_v <- path.expand(file.path(outDir_v, paste0("allClones_", wkbkName_v)))
summaryFile_v <- path.expand(file.path(outDir_v, paste0("summary_", wkbkName_v)))
plotDir_v <- mkdir(outDir_v, "pairwisePlots")
#plotFile_v <- path.expand(file.path(outDir_v, plotWkbkName_v))

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
  subMeta_dt <- subMeta_dt[get(compare_v) %in% whichCompare_v,]
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
                                currMeta_dt[[compare_v]], currMeta_dt[[sepCol_v]], sep = "_")
    
    ## Subset data
    sepData_lsdt[[val_v]] <- inputData_dt[Sample %in% currSamples_v,]
    sepMeta_lsdt[[val_v]] <- currMeta_dt
    
  } # for val_v
} # fi

####################
### PREPARE PLOT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

### Get divisions and colors
divisions_v <- c("Rare" = 0.00001, "Small" = 0.0001, "Medium" = 0.001, "Large" = 0.01, "Hyperexpanded" = 1)
divColors_v <- rev(brewer.pal(11, "RdYlBu")[c(1,3,5,9,11)])
names(divColors_v) <- names(divisions_v)

### Create data for rects
rects <- data.frame(ystart = c(0, divisions_v[1:4]), yend = divisions_v, col= names(divisions_v))
rects$col <- factor(rects$col, levels = rev(c("Rare", "Small", "Medium", "Large", "Hyperexpanded")))

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
  
  ## Make summary tables
  currFullSummary_dt <- currFullMerge_dt <- NULL
  
  ## Run for each animal
  for (j in 1:length(currMice_v)) {
    
    ###
    ### GET SPECIFIC MOUSE INFO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###
    
    ## Get mouse
    currMouse_v <- currMice_v[j]
    
    ## Subset metadata and get samples
    currMouseMeta_dt <- currMeta_dt[get(pairCol_v) == currMouse_v,]
    
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
      currMouseCompare_v <- currMouseMeta_dt[k,get(compare_v)]
      
      ## Make name including mouse, sample, and division 2
      currMouseName_v <- paste(currMouse_v, currMouseSample_v, currMouseCompare_v, sep = "_")
      
      ## If multiple sepData (i.e. using division 1) add that
      if (currName_v != "all") currMouseName_v <- paste(currMouseName_v, currName_v, sep = "_")
      
      ## Check that it matches 'master'
      if (currMouseName_v != currMouseMeta_dt[k,master]) stop("Mouse Name does not match 'master' name. Please check.")
      
      ## Subset for frequency groups and add to list
      currMouseData_lsdt[[currMouseName_v]] <- currData_dt[Sample == currMouseSample_v &
                                                             get(compare_v) == currMouseCompare_v,]
    } # for k
    
    ###
    ### SUBSET FOR SHARED CLONES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###
    
    ### Get all pairwise comparisons
    currMouseCompareVals_v <- sapply(names(currMouseData_lsdt), function(x) strsplit(x, split = "_")[[1]][3])
    currMouseCompareVals_v <- sort(currMouseCompareVals_v)
    currMouseCombn_lsv <- combn(currMouseCompareVals_v, 2, simplify = F)
    
    ### Get pairwise intersections
    currIntersect_lsv <- lapply(currMouseCombn_lsv, function(x) {
      intersect(currMouseData_lsdt[[ names(x)[1] ]]$compare, currMouseData_lsdt[[ names(x)[2] ]]$compare)
    })
    names(currIntersect_lsv) <- sapply(currMouseCombn_lsv, function(x) paste(x, collapse = "_"))
    
    ### Get total intersections
    currTotalIntersect_v <- unique(unlist(currIntersect_lsv))
    
    ### skip, if empty
    if (length(currTotalIntersect_v) == 0) {
      cat(sprintf("No overlap for mouse: %s\n", currMouse_v))
      next
    }
    
    ### Subset for clones
    currSharedClones_lsdt <- lapply(currMouseData_lsdt, function(x) {
      y <- x[x$compare %in% currTotalIntersect_v,]
      return(y)
    })
    
    ##
    ## OUTPUT TABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##
    
    ## Combine data 
    currSharedMerge_dt <- mergeDTs(currSharedClones_lsdt, mergeCol_v = "compare", keepCol_v = c("Div"))
    
    ## Change NA to ""
    for (col_v in names(currSharedClones_lsdt)) currSharedMerge_dt[is.na(get(col_v)), (col_v) := ""]
    
    ## Summarized
    currSharedSummary_dt <- as.data.table(table(currSharedMerge_dt[,mget(names(currSharedClones_lsdt))]))
    
    ## Add frequencies
    currSharedMerge_dt <- suppressWarnings(mergeDTs(currSharedClones_lsdt, mergeCol_v = "compare", keepCol_v = c("Div", "nb.clone.fraction")))
    
    ## Add sample columns to merge
    currSharedMergeDivCols_v <- grep("_Div", colnames(currSharedMerge_dt), value = T)
    currSharedMergeFreqCols_v <- grep("_nb.clone.fraction", colnames(currSharedMerge_dt), value = T)
    for (m in 1:length(currSharedMergeDivCols_v)) {
      newCol_v <- paste0('mouse', m)
      currSharedMerge_dt[, (newCol_v) := rep(currSharedMergeDivCols_v[m], nrow(currSharedMerge_dt))]
      colnames(currSharedMerge_dt)[grep(currSharedMergeDivCols_v[m], colnames(currSharedMerge_dt))] <- paste0("mouse", m, "_Div")
      colnames(currSharedMerge_dt)[grep(currSharedMergeFreqCols_v[m], colnames(currSharedMerge_dt))] <- paste0("mouse", m, "_Freq")
    }
    
    ## Add sample columns to summary
    for (m in 1:length(names(currSharedClones_lsdt))) {
      newCol_v <- paste0('mouse', m)
      currSharedSummary_dt[, (newCol_v) := rep(names(currSharedClones_lsdt)[m], nrow(currSharedSummary_dt))]
      colnames(currSharedSummary_dt)[colnames(currSharedSummary_dt) == names(currSharedClones_lsdt)[m]] <- paste0("mouse", m, "_Div")
    }
    
    ## Add third mouse, if needed
    currSharedOutMerge_dt <- currSharedMerge_dt; currSharedOutMerge_dt[1,1] <- currSharedOutMerge_dt[1,1]
    currSharedOutSummary_dt <- currSharedSummary_dt; currSharedOutSummary_dt[1,1] <- currSharedOutSummary_dt[1,1]
    if (length(currSharedClones_lsdt) < 3) {
      currSharedOutMerge_dt$mouse3_Div <- currSharedOutMerge_dt$mouse3_Freq <- currSharedOutMerge_dt$mouse3 <- rep(NA, nrow(currSharedOutMerge_dt))
      currSharedOutSummary_dt$mouse3_Div <- currSharedOutSummary_dt$mouse3 <- rep(NA, nrow(currSharedOutSummary_dt))
    }
    
    ## Fix column order
    mergeOutCols_v <- c("compare", 
                        grep("^mouse[1-9]$", colnames(currSharedOutMerge_dt), value = T),
                        grep("_Div$", colnames(currSharedOutMerge_dt), value = T), 
                        grep("_Freq$", colnames(currSharedOutMerge_dt), value = T))
    summaryOutCols_v <- c(grep("^mouse[1-9]$", colnames(currSharedOutSummary_dt), value = T),
                          grep("_Div$", colnames(currSharedOutSummary_dt), value = T), "N")
    
    ## Add to summary
    currFullMerge_dt <- rbind(currFullMerge_dt, currSharedOutMerge_dt[,mget(mergeOutCols_v)])
    currFullSummary_dt <- rbind(currFullSummary_dt, currSharedOutSummary_dt[,mget(summaryOutCols_v)])
    
    ##
    ## OUTPUT GRAPH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##
    
    ## Melt
    idCols_v <- grep("_Freq", colnames(currSharedMerge_dt), value = T, invert = T)
    currSharedMelt_dt <- melt(currSharedMerge_dt, id.vars = idCols_v)
    
    ## Title and file
    currTitle_v <- paste0("Mouse - ", currMouse_v)
    currFile_v <- paste0("mouse_", currMouse_v, "_trackClones.pdf")
    if (currName_v != "all") {
      currTitle_v <- paste0(currTitle_v, " (", currName_v, ")")
      currFile_v <- paste0(currName_v, "_", currFile_v)
    }
    
    ## Get xmax and xmin
    all_v <- unique(as.character(currSharedMelt_dt$variable))
    xmin_v <- all_v[which.min(as.numeric(gsub("mouse|_Freq", "", all_v)))]
    xmax_v <- all_v[which.max(as.numeric(gsub("mouse|_Freq", "", all_v)))]
    
    ## Add names for labels
    names(all_v) <- currMouseCompareVals_v
    
    ## Add color variable
    currSharedMelt_dt$Color <- apply(currSharedMelt_dt, 1, function(x) {
      if (length(which(is.na(x))) > 0) {return("sub")} else {return("all")}
    })
    
    ## Remove NA values (this is so that if a clone is in time1 and time3, but not time2, it will still get a line)
    currSharedMelt_dt <- currSharedMelt_dt[!is.na(value),]
    
    
    ## Plot - with color
    currPlot_gg <- ggplot() +
      geom_rect(data = rects, aes(xmin = xmin_v, xmax = xmax_v, ymin = ystart, ymax = yend, fill = col), alpha = 0.5) +
      scale_fill_manual(values = divColors_v, labels = rev(names(divColors_v))) +
      # geom_line(data = currSharedMelt_dt, aes(x = variable, y = value, group = compare, color = Color), alpha = 0.4) +
      geom_line(data = currSharedMelt_dt, aes(x = variable, y = value, group = compare, color = Color)) +
      geom_point(data = currSharedMelt_dt, aes(x = variable, y = value, group = compare, color = Color), alpha = 0.4) +
      scale_color_manual(values = c(alpha("black", 0.4), alpha("darkgrey", 0.6)), labels = c("all", "sub")) +
      scale_y_log10(breaks = unname(divisions_v)) +
      scale_x_discrete(breaks = all_v, labels = names(all_v))+
      ggtitle(currTitle_v) +
      big_label() +
      labs(y = "Normalized Clone Frequency", fill = "Freq Group") +
      guides(color = F) +
      theme(axis.title.x = element_blank())
    
    ## Save plot
    ggsave(filename = file.path(plotDir_v, currFile_v), plot = currPlot_gg, width = 7, height = 7)
    
      
  } # for j
  
  ##
  ## WRITE OUTPUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##
  
  ## Write - this needs to be changed to fit.
  if ("all" %in% toWrite_v) {
    append_v <- file.exists(allFile_v)
    write.xlsx(x = currFullMerge_dt, file = allFile_v, sheetName = paste(currName_v, sheetName_v, "-_-"),
               row.names = F, append = append_v)
  } # fi

  if ("summary" %in% toWrite_v) {
    append_v <- file.exists(allFile_v)
    write.xlsx(x = currFullSummary_dt, file = summaryFile_v, sheetName = paste(currName_v, sheetName_v, "-_-"),
               row.names = F, append = append_v)
  }
  
} # for i

