#!/usr/bin/Rscript

###################
### DESCRIPTION ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################

### Generate overall read count and contamination QC plots
  ### Plot 1
    ### Violin/boxplot of total reads for all samples.
    ### Add points and label points with < 100k reads
  ### Plot 2
    ### Same as above, but for clones.
  ### Plot 3/4
    ### Compare total reads to total clones and uniq clones
  ### Plot 5
    ### "Scatter" plot of rank/number of clones for contaminants
    ### X-axis is sample
    ### Y-axis is rank/numClones
    ### Each sample has 3 points, one for each contaminant
  ### Plot 6
    ### Boxplot of rank/numClones
    ### X-axis is contaminant
    ### Y-axis is rank/numClones
    ### Each box (there will be 3) has 1 data point for each sample
  ### Plot 7
    ### Percentage of each sample occupied by contaminants (in terms of clones)
  ### Contamination Summary
    ### Record how many (and which) samples have rank <= 15 for at least one contaminant
    ### Record how many (and which) samples have rank <= 15, summarized by contaminant.

####################
### DEPENDENCIES ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(data.table))
suppressMessages(library(optparse))
library(wrh.rUtils)
options("scipen" = 100)
options(warn = 0)

####################
### COMMAND LINE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

### Make command line
optlist <- list(
  make_option(
    c("-C", "--contamFile"),
    type = "character",
    help = "Path to contaminationQC file. Usually called: LIB######LC_contaminationQC.txt"
  ),
  make_option(
    c("-N", "--nineFile"),
    type = "character",
    help = "Path to 9bp count spikes QC summary table. Should be: count.spikes.9bp.QC.summary.txt"
  ),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    help = "Path to directory to write output files. If this argument is not used, will write to location of inFile."
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
    Default is 'Treatment'. Other common columns could be Timepoint.
    If metaFile is provided and you would not like to separate by color, provide 'NA' to skip."
  ),
  make_option(
    c("-f", "--facetCol"),
    type = "character",
    help = "(Optional) metadata column used to separate samples by facets. Only used if metaFile is provided.
    Default is 'Tissue'. Other common columns could be Type (mouse type), or Tumor (tumor/non-tumor)
    If metaFile is provided and you would not like to facet, provide 'NA' to skip."
  )
)

### Parse command line
p <- OptionParser(usage = "%prog -C contamFile -N nineFile -o outDir -m metaFile -c colorCol -f facetCol",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

### Get arguments
contamFile_v <- args$contamFile
nineFile_v <- args$nineFile
outDir_v <- args$outDir
metaFile_v <- args$metaFile
colorCol_v <- args$colorCol
facetCol_v <- args$facetCol

### For testing
# contamFile_v <- "~/OHSU/tcr_spike/data/LIB190308LC/qc/LIB190308LC_contaminationQC.txt"
# nineFile_v <- "~/OHSU/tcr_spike/data/LIB190308LC/qc/count.spikes.9bp.QC.summary.txt"
# outDir_v <- mkdir("~/OHSU/tcr_spike/data/LIB190308LC/qc_plots/", "reads_and_contam")
# metaFile_v <- NULL
# colorCol_v <- NULL
# facetCol_v <- NULL

###################
### HANDLE NULL ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################

### outDir_v becomes inputDir_v, if NULL
if (is.null(outDir_v)) {
  outDir_v <- dirname(contamFile_v)
  cat(sprintf("'--outDir' not specified. Writing output to location of contamFile:\n %s", outDir_v))
}

### Lots of metadata stuff to do if it's provided. Do nothing otherwise
if (!is.null(metaFile_v)) {
  
  ## Read data
  meta_dt <- fread(metaFile_v)
  
  ## Get sample column
  sampleCol_v <- grep("[Ss]ample", colnames(meta_dt), value = T)
  if (length(sampleCol_v) == 0) {
    stop("metaFile provided, but couldn't find 'Sample' column.\n",
         sprintf("Current columns are: %s\n", paste(colnames(meta_dt), collapse = " - ")),
         "Please rename the column containing sample numbers to 'Sample'.")
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

#################
### PREP DATA ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################

### Read data
contam_dt <- fread(contamFile_v)
nine_dt <- fread(nineFile_v)

### Fix spike count data
nine_dt <- nine_dt[,mget(grep("DM_*", colnames(nine_dt), value = T, invert = T))]
nine_dt$sample.id <- paste0("S", gsub("^.*_S|\\..*$", "", nine_dt$sample.id))

### Add read count flag (flag all reads w/ fewer than 10% of mean)
medianReads_v <- median(nine_dt$num.reads)
medianCut_v <- medianReads_v * 0.1
nine_dt$flag <- ""
nine_dt[num.reads < medianCut_v, flag := sample.id]

### Add total and unique clone flags
medianTotalClone_v <- median(contam_dt$Orig.Total.Clones); medianTotalCut_v <- medianTotalClone_v * 0.1
medianUniqClone_v <- median(contam_dt$Orig.Unique.Clones); medianUniqCut_v <- medianUniqClone_v * 0.1
contam_dt$totalFlag <- contam_dt$uniqFlag <- ""
contam_dt[Orig.Total.Clones < medianTotalCut_v, totalFlag := Sample]
contam_dt[Orig.Unique.Clones < medianUniqCut_v, uniqFlag := Sample]

### Add dummy variables
nine_dt$dummy <- "1"
contam_dt$dummy <- "1"

########################
### READ/CLONE PLOTS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################

### Total Reads
totalReadRaw_gg <- ggplot(nine_dt, aes(y = num.reads, x = dummy)) + 
  geom_violin() +
  geom_point() + 
  ggrepel::geom_text_repel(aes(y = num.reads, label = flag)) +
  labs(x = NULL, y = "Total Reads") +
  ggtitle("Distribution of Reads\n<10% median are labeled") +
  big_label() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

totalReadLog_gg <- totalReadRaw_gg + scale_y_log10()

### Total Clones
totalCloneRaw_gg <- ggplot(contam_dt, aes(y = Orig.Total.Clones, x = dummy)) +
  geom_violin() +
  geom_point() +
  ggrepel::geom_text_repel(aes(y = Orig.Total.Clones, label = totalFlag)) +
  labs(x = NULL, y = "Total Clones") +
  ggtitle("Distribution of Total Clones\n<10% median are labeled") +
  big_label() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
totalCloneLog_gg <- totalCloneRaw_gg + scale_y_log10()

### Unique Clones
uniqCloneRaw_gg <- ggplot(contam_dt, aes(y = Orig.Unique.Clones, x = dummy)) +
  geom_violin() +
  geom_point() +
  ggrepel::geom_text_repel(aes(y = Orig.Unique.Clones, label = totalFlag)) +
  labs(x = NULL, y = "Unique Clones") +
  ggtitle("Distribution of Unique Clones\n<10% median are labeled") +
  big_label() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
uniqCloneLog_gg <- uniqCloneRaw_gg + scale_y_log10()

### Print
plots_lslsgg <- list("raw" = list("totalReads" = totalReadRaw_gg, "totalClones" = totalCloneRaw_gg, "uniqClones" = uniqCloneRaw_gg),
                     "log10" = list("totalReads" = totalReadLog_gg, "totalClones" = totalCloneLog_gg, "uniqClones" = uniqCloneLog_gg))

for (i in 1:length(plots_lslsgg)) {
  name1_v <- names(plots_lslsgg)[[i]]
  for (j in 1:length(plots_lslsgg[[name1_v]])) {
    name2_v <- names(plots_lslsgg[[name1_v]])[j]
    ggsave(filename = file.path(outDir_v, paste(name1_v, name2_v, "distr.pdf", sep = "_")), 
           plot = plots_lslsgg[[name1_v]][[name2_v]], width = 5, height = 5)
  }
}

##############################
### READ/CLONE CORRELATION ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############################

### Merge data
merge_dt <- merge(nine_dt[,mget(c("sample.id", "num.reads"))], 
                  contam_dt[,mget(c("Sample", "Orig.Unique.Clones", "Orig.Total.Clones"))],
                  by.x = "sample.id", by.y = "Sample", sort = F)

### Label
merge_dt$lab <- ""
#merge_dt[num.reads < 250000, lab := sample.id]
merge_dt[(num.reads < 50000 & Orig.Total.Clones < 100000), lab := sample.id]

### Variables
vars_v <- c("Orig.Unique.Clones", "Orig.Total.Clones")

### Plot
for (i in 1:length(vars_v)) {
  
  ## Get var and label
  currVar_v <- vars_v[i]
  currVarLab_v <- ifelse(currVar_v == "Orig.Unique.Clones", "Unique Clones", "Total Clones")
  
  ## Get formula
  currLM <- lm(merge_dt[[currVar_v]] ~ merge_dt$num.reads)
  currR2 <- round(summary(currLM)$adj.r.squared, digits = 3)
  currP <- modelP(currLM)
  
  ## Make plot
  currPlot_gg <- ggplot(data = merge_dt, aes_string(x = "num.reads", y = currVar_v)) +
    geom_point() +
    ggrepel::geom_text_repel(aes_string(x = "num.reads", y = currVar_v, label = "lab")) +
    geom_smooth(method = 'lm', se = F, show.legend = T) +
    annotate("text", x = min(merge_dt$num.reads), y = max(merge_dt[[currVar_v]]), label = paste0("R2 = ", currR2), hjust = -1) +
    labs(x = "Total Reads", y = currVarLab_v) +
    ggtitle(paste0("Compare Reads to ", currVarLab_v)) +
    scale_y_log10() + scale_x_log10() +
    big_label()
  
  ## Print
  ggsave(filename = file.path(outDir_v, paste0("correlationReadsTo_", gsub(" ", "_", currVarLab_v), ".pdf")),
         plot = currPlot_gg, width = 5, height = 5)
  
}

#####################################
### CONTAMINATION PLOTS - SCATTER ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####################################

### Melt
contamMelt_dt <- melt(contam_dt, id.vars = c("Sample", "uniqFlag", "totalFlag", "dummy"))

### Variables
vars_lsv <- list("Rank" = c("p14.rank", "ot1.rank", "el4.rank"),
                 "Count" = c("p14.count", "ot1.count", "el4.count"))

### Plot
for (i in 1:length(vars_lsv)) {
  
  ###
  ### PREP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ## Get subset info
  currVarType_v <- names(vars_lsv)[i]
  currVars_v <- vars_lsv[[currVarType_v]]
  
  ## Subset
  currData_dt <- contamMelt_dt[variable %in% currVars_v,]
  
  ## Fix variable names (i.e. remove 'rank' or 'count)
  currData_dt$variable <- gsub("\\.rank|\\.count", '', currData_dt$variable)
  
  ## Sort samples
  currData_dt <- currData_dt[order(value),]
  currData_dt$Sample <- factor(currData_dt$Sample, levels = unique(currData_dt$Sample))
  
  ###
  ### ALL SAMPLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ## Plot
  currPlot_gg <- ggplot(currData_dt, aes(x = Sample, y = value, color = variable)) +
    geom_point() +
    scale_y_log10() +
    labs(y = currVarType_v, color = "Contaminant") +
    ggtitle(paste0(currVarType_v, " of Contaminant Sequences")) +
    big_label() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  ###
  ### CLOSE UP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ## Subset
  if (currVarType_v == "Rank") {
    close_v <- 15
    currClose_dt <- currData_dt[value < close_v,]
  } else {
    close_v <- quantile(currData_dt[value != 0,value], probs = 0.9)
    currClose_dt <- currData_dt[value > close_v,]
  }
  
  ## Text size
  textSize_v <- ifelse(nrow(currClose_dt) > 20, 7, 15)
  
  ## Plot
  currClose_gg <- ggplot(currClose_dt, aes(x = Sample, y = value, color = variable)) +
    geom_point() +
    scale_y_log10() +
    labs(y = currVarType_v, color = "Contaminant") +
    ggtitle(paste0(currVarType_v, " of Contaminant Sequences")) +
    big_label() + angle_x() +
    theme(axis.text.x = element_text(size = textSize_v))
  currClose_gg
  
  ###
  ### PRINT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Full
  ggsave(filename = file.path(outDir_v, paste0("contam_", currVarType_v, ".pdf")),
         plot = currPlot_gg, width = 7, height = 7)
  
  ### Close up
  ggsave(filename = file.path(outDir_v, paste0("topContam_", currVarType_v, ".pdf")),
         plot = currClose_gg, width = 7, height = 7)
  
} # for i

#################################
### CONTAMINATION PLOTS - BOX ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################################

### Plot
for (i in 1:length(vars_lsv)) {
  
  ###
  ### PREP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ## Get subset info
  currVarType_v <- names(vars_lsv)[i]
  currVars_v <- vars_lsv[[currVarType_v]]
  
  ## Subset
  currData_dt <- contamMelt_dt[variable %in% currVars_v,]
  
  ## Fix variable names (i.e. remove 'rank' or 'count)
  currData_dt$variable <- gsub("\\.rank|\\.count", '', currData_dt$variable)
  
  ## Add labels
  currData_dt$lab <- ""
  if (currVarType_v == "Rank") {
    close_v <- 10
    currData_dt[value < close_v, lab := Sample]
  } else {
    close_v <- quantile(currData_dt[value != 0,value], probs = 0.9)
    currData_dt[value > close_v, lab := Sample]
  }
  
  ###
  ### PLOT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###

  currPlot_gg <- ggplot(currData_dt, aes(x = variable, y = value, color = variable)) +
    geom_violin() +
    geom_point() +
    scale_y_log10() +
    ggrepel::geom_text_repel(aes(y = value, label = lab)) +
    labs(y = currVarType_v, x = "Contaminant") +
    guides(color = F) +
    ggtitle(paste0("Distr. of Contaminant ", currVarType_v)) +
    big_label()
  
  ###
  ### PRINT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Full
  ggsave(filename = file.path(outDir_v, paste0("contaminantDistr_", currVarType_v, ".pdf")),
         plot = currPlot_gg)
  
} # for i

#####################################
### CONTAMINATION PLOTS - PERCENT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####################################

### Calculate percentages
cols_v <- c("all" = "Contam.clones", "p14" = "p14.count", "ot1" = "ot1.count", "el4" = "el4.count")

for (i in 1:length(cols_v)) {
  
  ## Get columns
  col_v <- cols_v[i]
  newCol_v <- names(cols_v)[i]
  
  ## Calculate
  contam_dt[[newCol_v]] <- contam_dt[[col_v]] / contam_dt$Orig.Total.Clones * 100
}

### Subset and melt
contamMelt_dt <- melt(contam_dt[,mget(c("Sample", names(cols_v)))], id.vars = "Sample")
contamMelt_dt$dummy <- "1"

### Add label
contamMelt_dt$lab <- ""
contamMelt_dt[value > 5, lab := Sample]

### Make Plot
pct_gg <- ggplot(data = contamMelt_dt, aes(x = variable, y = value, color = variable)) +
  geom_violin() +
  geom_point() +
  scale_y_log10(breaks = c(0.01, 1, 10, 100)) +
  guides(color = F) +
  ggrepel::geom_text_repel(aes(y = value, label = lab)) +
  ggtitle("Percentage of Total Clones in Contaminants") +
  labs(y = "Percent of All Clones", x = "Contaminant") +
  big_label()

### Print
ggsave(filename = file.path(outDir_v, "percentClonesContaminated.pdf"), plot = pct_gg)