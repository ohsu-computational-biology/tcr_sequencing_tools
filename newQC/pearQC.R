#!/usr/bin/Rscript

###################
### DESCRIPTION ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################

### 03_pear.pl sometimes has trouble producing the pear_summary_log file.
### Each sample produces a file called "pear_full_log_[#].txt that contains
### detailed output of the run.
###
### Read each of these files and create a new pear_summary_log.txt file
### with one row for each sample and 3 data columns
    ### pct Assembled
    ### pct Discarded
    ### pct Unassembled
###
### Also create a violin (or box) plot for each of the 3 columns.
### If a metadata file is provided, split the violins by color and/or facet.
### Color generally corresponds to treatment, while facet corresponds to tissue, but can change.

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
    c("-i", "--inputDir"),
    type = "character",
    help = "Path to directory of pear log files. Will only read files that satisfy pear_full_log_[0-9]+.txt"
  ),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    help = "Path to directory to write output files. If this argument is not used, will write to inputDir."
  ),
  make_option(
    c("-m", "--metaFile"),
    type = "character",
    help = "(Optional) path to metadata file. Used to divide plots by treatment, if provided."
  ),
  make_option(
    c("-c", "--colorCol"),
    type = "character",
    help = "(Optional) metadata column name used to separate samples using different colors. Only used if metaFile is provided.
    Default is 'Treatment'. Provide 'NA' to skip."
  ),
  make_option(
    c("-f", "--facetCol"),
    type = "character",
    help = "(Optional) metadata column used to separate samples by facets. Only used if metaFile is provided.
    Default is 'Tissue'. Provide 'NA' to skip."
  )
)

### Parse command line
p <- OptionParser(usage = "%prog -i inputDir -o outDir -m metaFile -c colorCol -f facetCol",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

### Get arguments
inputDir_v <- args$inputDir
outDir_v <- args$outDir
metaFile_v <- args$metaFile
colorCol_v <- args$colorCol
facetCol_v <- args$facetCol

### For testing
# inputDir_v <- "~/OHSU/tcr_spike/data/test_collab/pear/"
# outDir_v <- "~/OHSU/tcr_spike/data/test_collab/out/"
# metaFile_v <- "~/OHSU/tcr_spike/data/test_collab/meta.txt"
# colorCol_v <- "Treatment"
# facetCol_v <- "Source"

###################
### HANDLE NULL ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################

### outDir_v becomes inputDir_v, if NULL
if (is.null(outDir_v)) outDir_v <- inputDir_v

### Lots of metadata stuff to do if it's provided. Do nothing otherwise
if (!is.null(metaFile_v)) {
  
  ## Read data
  meta_dt <- fread(metaFile_v)
  
  ## Get sample column
  sampleCol_v <- grep("[Ss]ample", colnames(meta_dt), value = T)
  if (length(sampleCol_v) == 0) {
    stop("metaFile provided, but couldn't find 'Sample' column.\n",
         sprintf("Current columns are: %s\n", paste(colnames(meta_dt), collapse = " - ")),
         "Please rename the column with sample numbers to 'Sample'.")
  } # fi
  
  ## Get color-separating column
  colorCol_v <- ifelse(is.null(colorCol_v), grep("[Tt]reatment|[Tt]reat", colnames(meta_dt), value = T), colorCol_v)
  if (is.na(colorCol_v)) {
    colorCol_v <- NULL
    warning("metaFile provided, but color column was NULL and no corresponding column was found.\n",
            sprintf("Current columns are: %s\n", paste(colnames(meta_dt), collapse = " - ")),
            "This is only a problem if you would like to split plots by a particular variable.\n")
  } else {
    cat(sprintf("Will split samples by color using: %s\n", colorCol_v))
  } # fi
  
  ## Get faceting column
  facetCol_v <- ifelse(is.null(facetCol_v), grep("[Tt]issue", colnames(meta_dt), value = T), facetCol_v)
  if (is.na(facetCol_v)) {
    facetCol_v <- NULL
    message("metaFile provided, but facet column was NULL and no corresponding column was found.\n",
            sprintf("Current columns are: %s\n", paste(colnames(meta_dt), collapse = " - ")),
            "This is only a problem if you would like to facet plots by a particular variable.\n")
  } else {
    cat(sprintf("Will facet plots by: %s\n", facetCol_v))
  } # fi
  
  ## Add 'S' to front of samples
  meta_dt[[sampleCol_v]] <- paste0("S", gsub("^S", "", meta_dt[[sampleCol_v]]))
  
} # fi

#############
### SETUP ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############

### Get files, add names, and sort them
inputFiles_v <- list.files(inputDir_v, pattern = "pear_full_log_[0-9]+.txt")
names(inputFiles_v) <- paste0("S", gsub("^.*_|\\.txt", "", inputFiles_v))
inputFiles_v <- inputFiles_v[order(as.numeric(gsub("S", "", names(inputFiles_v))))]

### Make output data.frame
cols_v <- paste0("pear_pct", c("Assembled", "Unassembled", "Discarded"))
out_df <- data.frame("Sample" = names(inputFiles_v), stringsAsFactors = F)
for (col_v in cols_v) out_df[[col_v]] <- numeric(length(inputFiles_v))

### Make output names
outFile_v <- file.path(outDir_v, "pear_summary_log.txt")
plotFiles_v <- sapply(cols_v, function(x) file.path(outDir_v, paste0(x, ".pdf")))

############
### READ ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############

for (i in 1:length(inputFiles_v)) {
  
  ## Get name and record
  currSample_v <- names(inputFiles_v)[i]
  currRecord <- readLines(file.path(inputDir_v, inputFiles_v[i]))
  
  ## Get batch name
  batch_v <- gsub("_.*$", "", basename(strsplit(grep("Forward reads file.*", currRecord, value = T), split = ":")[[1]][2]))
  
  ## Get percentages
  currAssemble_v <- as.numeric(strsplit(grep("Assembled reads \\.", currRecord, value = T), split = "\\(|%")[[1]][2])
  currUnassemble_v <- as.numeric(strsplit(grep("Not assembled reads \\.", currRecord, value = T), split = "\\(|%")[[1]][2])
  currDiscarded_v <- as.numeric(strsplit(grep("Discarded reads \\.", currRecord, value = T), split = "\\(|%")[[1]][2])
  
  ## Add to output
  out_df$Sample[i] <- currSample_v
  out_df[i,cols_v] <- c(currAssemble_v, currUnassemble_v, currDiscarded_v)
}

write.table(out_df, outFile_v, sep = '\t', quote = F, row.names = F)

##############
### NOTIFY ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############

for (i in 1:nrow(out_df)) {
  currAssemble_v <- out_df[i,"pear_pctAssembled"]
  if (currAssemble_v < 90) {
    cat(sprintf("Sample %s only assembled %s%% of reads.\n", out_df[i,"Sample"], as.character(currAssemble_v)))
  } # fi
} # for

############
### PLOT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############

if (is.null(metaFile_v)) {
  
  ## Melt
  melt_df <- melt(out_df, id.vars = "Sample")
  
  for (i in 1:length(cols_v)) {
    
    ## Get info
    currCol_v <- cols_v[i]
    currData_df <- melt_df[melt_df$variable == currCol_v, ]
    currOutFile_v <- plotFiles_v[i]
    
    ## Make plot
    currPlot_gg <- ggplot(data = currData_df, aes(x = variable, y = value)) +
      geom_violin(trim = F) +
      stat_summary(fun.y = median, geom = "point", size = 2) +
      ylab("Percent") +
      ggtitle(paste(batch_v, "Percent Reads", gsub("pct", "", currCol_v))) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 20),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 16),
            axis.title.y = element_text(size = 18))
    
    ## Write plot
    pdf(file = currOutFile_v)
    print(currPlot_gg)
    dev.off()
  } # for i
} else {
  
  ## Merge 
  metaCols_v <- c(colorCol_v, facetCol_v)
  merge_dt <- merge(meta_dt[,mget(c(sampleCol_v, metaCols_v))], out_df, by.x = sampleCol_v, by.y = "Sample", sort = F)
  
  ## Melt
  melt_dt <- melt(merge_dt, id.vars = c(sampleCol_v, metaCols_v))
  
  for (i in 1:length(cols_v)) {
    
    ## Get info
    currCol_v <- cols_v[i]
    currData_dt <- melt_dt[melt_dt$variable == currCol_v, ]
    currOutFile_v <- plotFiles_v[i]
    
    ## Determine if violin or box plot
    ## Box plot if at least one grouping has only 2 observations
    sub_dt <- currData_dt[,mget(metaCols_v)]
    count_table <- table(sub_dt)
    minObs_v <- min(count_table)
    
    ## Make base plot
    currPlot_gg <- ggplot(data = currData_dt, aes(x = variable, y = value))
    
    ## Add either boxplot or violin plot
    if (minObs_v > 2) {
      currPlot_gg <- currPlot_gg + geom_violin()
    } else {
      currPlot_gg <- currPlot_gg + geom_boxplot()
    } # fi
    
    ## Split by color
    if (!is.null(colorCol_v)) {
      currPlot_gg <- currPlot_gg + aes_string(color = colorCol_v)
    }
    
    ## Facet
    if (!is.null(facetCol_v)) {
      currPlot_gg <- currPlot_gg + facet_wrap(as.formula(paste("~", facetCol_v)))
    }
    
    ## Add other info
    currPlot_gg <- currPlot_gg +
      ylab("Percent") +
      ggtitle(paste(batch_v, "Percent Reads", gsub("pct", "", currCol_v)))
    
    ## Modify theme
    currPlot_gg <- currPlot_gg +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 20),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 16),
            axis.title.y = element_text(size = 18),
            strip.text = element_text(size = 18),
            legend.position = "bottom")

    ## Write plot
    pdf(file = currOutFile_v)
    print(currPlot_gg)
    dev.off()
  } # for i

} # fi
