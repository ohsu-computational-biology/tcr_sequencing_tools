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
suppressMessages(library(ggpubr))
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
  ),
  make_option(
    c("-p", "--plot"),
    type = "character",
    help = "'dot' = dot plot; 'box' = box plot"
  ),
  make_option(
    c("-t", "--test"),
    type = "logical",
    default = F,
    help = "TRUE - run statistical test; FALSE - no test."
  )
)

### Parse command line
p <- OptionParser(usage = "%prog -C contamFile -N nineFile -o outDir -m metaFile -c colorCol -f facetCol -p plot -t test",
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
plot_v <- args$plot
compare_v <- args$test

### For testing
# contamFile_v <- "~/OHSU/tcr_spike/data/LIB190701LC/compareShortLong/input/LIB190701LC_contaminationQC.txt"
# nineFile_v <- "~/OHSU/tcr_spike/data/LIB190701LC/compareShortLong/input/count.spikes.9bp.QC.summary.txt"
# outDir_v <- mkdir("~/OHSU/tcr_spike/data/LIB190701LC/compareShortLong/", "output")
# metaFile_v <- "~/OHSU/tcr_spike/data/LIB190701LC/meta/meta.txt"
# colorCol_v <- "Experiment"
# facetCol_v <- "Tissue"
# plot_v <- "box"
# compare_v <- T

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

### Merge with meta
if (!is.null(metaFile_v)) {
  nine_dt <- merge(nine_dt, meta_dt, by.x = "sample.id", by.y = "Sample", sort = F)
  contam_dt <- merge(contam_dt, meta_dt, by = "Sample", sort = F)
}

########################
### READ/CLONE PLOTS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################

### BASE PLOTS
totalReadRaw_gg <- ggplot(nine_dt, aes(y = num.reads, x = dummy)) +
  labs(x = NULL, y = "Total Reads") +
  ggtitle("Distr. of Reads") +
  big_label() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

totalCloneRaw_gg <- ggplot(contam_dt, aes(y = Orig.Total.Clones, x = dummy)) +
  labs(x = NULL, y = "Total Clones") +
  ggtitle("Distr. of Total Clones") +
  big_label() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

uniqCloneRaw_gg <- ggplot(contam_dt, aes(y = Orig.Unique.Clones, x = dummy)) +
  labs(x = NULL, y = "Unique Clones") +
  ggtitle("Distr. of Unique Clones") +
  big_label() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

### Select plot type
if (plot_v == "dot") {
  totalReadRaw_gg <- totalReadRaw_gg + geom_point(position = position_dodge(0.1))
  totalCloneRaw_gg <- totalCloneRaw_gg + geom_point(position = position_dodge(0.1))
  uniqCloneRaw_gg <- uniqCloneRaw_gg + geom_point(position = position_dodge(0.1))
} else if (plot_v == 'box') {
  totalReadRaw_gg <- totalReadRaw_gg + geom_boxplot()
  totalCloneRaw_gg <- totalCloneRaw_gg + geom_boxplot()
  uniqCloneRaw_gg <- uniqCloneRaw_gg + geom_boxplot()
}

### Add Colors
if (!is.null(colorCol_v)) {
  totalReadRaw_gg <- totalReadRaw_gg + aes_string(color = colorCol_v)
  totalCloneRaw_gg <- totalCloneRaw_gg + aes_string(color = colorCol_v)
  uniqCloneRaw_gg <- uniqCloneRaw_gg + aes_string(color = colorCol_v)
}

### Add facets
if (!is.null(facetCol_v)) {
  fxn <- paste0("~", facetCol_v)
  totalReadRaw_gg <- totalReadRaw_gg + facet_wrap(fxn)
  totalCloneRaw_gg <- totalCloneRaw_gg + facet_wrap(fxn)
  uniqCloneRaw_gg <- uniqCloneRaw_gg + facet_wrap(fxn)
}

### Compare
if (compare_v) {
  totalReadRaw_gg <- totalReadRaw_gg + stat_compare_means(aes(label = ..p.signif..))
  totalCloneRaw_gg <- totalCloneRaw_gg + stat_compare_means(aes(label = ..p.signif..))
  uniqCloneRaw_gg <- uniqCloneRaw_gg + stat_compare_means(aes(label = ..p.signif..))
}

### Add points
if (is.null(colorCol_v) & is.null(facetCol_v)) {
  totalReadRaw_gg <- totalReadRaw_gg + geom_point() + ggrepel::geom_text_repel(aes(y = num.reads, label = flag))
  totalCloneRaw_gg <- totalCloneRaw_gg + geom_point() + ggrepel::geom_text_repel(aes(y = Orig.Total.Clones, label = totalFlag))
  uniqCloneRaw_gg <- uniqCloneRaw_gg + geom_point() + ggrepel::geom_text_repel(aes(y = Orig.Unique.Clones, label = totalFlag))
}

### Make log10 versions
totalReadLog_gg <- totalReadRaw_gg + scale_y_log10()
totalCloneLog_gg <- totalCloneRaw_gg + scale_y_log10()
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

