#	RScript that calculates number of unique clonotypes, Shannon entropy, and clonality for multiple clonotype files
#
#   Input is a directory containing one or more files with MiXCR output format
#   
#   Load necessary libraries
library(entropy);
library(data.table);
library(tcR);
library(ggplot2);

#	Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);
#   Directory should contain all and only .txt files from exportClones, post-normalization with the following
#   format:
#     Clone count   Clone fraction    Clonal sequence(s)    AA. Seq. CDR3   Best V Hit    Best J Hit    V segments
#     J segments    Normalized clone count    Normalized clone fraction

cloneDir_v <- arguments[1];    # Typically .../dhaarini/DNAXXXXLC/normalization/normalized_clones/
metadata_v <- arguments[2]
outDir_v <- arguments[3];

### Get files and names
cloneFiles_v <- list.files(cloneDir_v)
cloneFiles_v <- cloneFiles_v[order(as.numeric(gsub("^.*_S|_.*$", "", cloneFiles_v)))]
cloneNames_v <- sapply(cloneFiles_v, function(x) strsplit(x, split = "_")[[1]][2], USE.NAMES=F)
batchName_v <- strsplit(cloneFiles_v[1], split = "_")[[1]][1]

### Read in data
cloneData_lsdt <- sapply(cloneFiles_v, function(x) fread(file.path(cloneDir_v, x)), simplify = F)
names(cloneData_lsdt) <- cloneNames_v

metadata_dt <- fread(metadata_v)

### Subset files

### For LIB170407LC
## metadata_dt <- metadata_dt[-c(4,21),]
## metadata_dt <- metadata_dt[tissue == "tumor",]
## keepSamp_v <- paste0("S", metadata_dt[, `sampleID/barcode`])
## cloneData_lsdt <- cloneData_lsdt[keepSamp_v]

### For LIB170213LC
## metadata_dt <- metadata_dt[c(1:35, 105:114),] #  graph 1
## metadata_dt <- metadata_dt[c(71:87, 125:134),] #  graph 2
## metadata_dt <- metadata_dt[c(88:104, 135:144),] #  graph 3

### Perform analysis on each treatment

## treatments_v <- unique(metadata_dt[, Treatment]) # LIB170407LC
## immuno_v <- unique(metadata_dt[, immuno]) # LIB170407LC
## treatments_v <- unique(metadata_dt[,treatment]) # LIB170213LC
treatments_v <- unique(metadata_dt[,Treatment]) # LIB170728LC

immunoHomeo_mat <- NULL
immunoRows_v <- NULL
finalHomeo_mat <- NULL
finalRows_v <- NULL

pdf(file = paste0(outDir_v, batchName_v, "_rmBadHomeostasis.pdf"))

for (i in 1:length(treatments_v)){

    ## Get Treatment
    currTreat_v <- treatments_v[i]

    ## Subset files
    ## currFiles_v <- paste0("S", metadata_dt[Treatment == currTreat_v, `sampleID/barcode`]) # LIB170407LC
    ## currFiles_v <- paste0("S", metadata_dt[treatment == currTreat_v, `sample`]) # LIB170213LC
    currFiles_v <- paste0("S", metadata_dt[Treatment == currTreat_v, `Sample`])
    currData_lsdt <- cloneData_lsdt[currFiles_v]

    ## Print nrow
    nrow_v <- sapply(currData_lsdt, nrow)
    print(c("Treatment" = currTreat_v, nrow_v))

    ## Remove zero-count clones
    currData_lsdt <- lapply(currData_lsdt, function(x) x[`Normalized clone fraction` > 0,])

    ## Run clonal space homeostasis
    currHomeo_mat <- clonal.space.homeostasis(.data = currData_lsdt, .prop.col = "Normalized clone fraction")

    ## Take mean of matrix
    currMean_v <- apply(currHomeo_mat, 2, mean)

    ## Add to list
    finalRows_v <- c(finalRows_v, currTreat_v)
    finalHomeo_mat <- rbind(finalHomeo_mat, currMean_v)

    ## Make plot
    currHomeo_gg <- vis.clonal.space(currHomeo_mat) +
        ggtitle(paste0(currTreat_v, " Clonal Space Homeostasis")) +
        theme(plot.title = element_text(hjust = 0.5))

    print(currHomeo_gg)
} # for i

rownames(finalHomeo_mat) <- finalRows_v

finalHomeo_gg <- vis.clonal.space(finalHomeo_mat) +
    ggtitle(paste0("Treat-wise Clonal Space Homeostasis")) +
    theme(plot.title = element_text(hjust = 0.5))

print(finalHomeo_gg)

### LIB170407LC ONLY
## null_mat <- finalHomeo_mat[1:2,]
## immuno_mat <- finalHomeo_mat[3:4,]

## final2Homeo_mat <- rbind(apply(null_mat, 2, mean), apply(immuno_mat, 2, mean))
## rownames(final2Homeo_mat) <- c("null", "aCD40")

## final2Homeo_gg <- vis.clonal.space(final2Homeo_mat) +
##     ggtitle("Immunotherapy Clonal Space Homeostasis") +
##     theme(plot.title = element_text(hjust = 0.5))

## print(final2Homeo_gg)

dev.off()

### LIB170407LC ONLY
## pdf(file = paste0(outDir_v, batchName_v, "_homeostasis2.pdf"))

## for (i in 1:length(immuno_v)){

##     ## Get immuno status
##     currImmuno_v <- immuno_v[i]

##     ## Subset files
##     currFiles_v <- paste0("S", metadata_dt[immuno == currImmuno_v, `sampleID/barcode`])
##     currData_lsdt <- cloneData_lsdt[currFiles_v]

##     ## Remove zero-count clones
##     currData_lsdt <- lapply(currData_lsdt, function(x) x[`Normalized clone fraction` > 0,])

##     ## Run clonal sapce homeostasis
##     currHomeo_mat <- clonal.space.homeostasis(.data = currData_lsdt, .prop.col = "Normalized clone fraction")

##     ## Take mean of matrix
##     currMean_v <- apply(currHomeo_mat, 2, mean)

##     ## Add to list
##     immunoRows_v <- c(immunoRows_v, currImmuno_v)
##     immunoHomeo_mat <- rbind(immunoHomeo_mat, currMean_v)

##     ## Make plot
##     currHomeo_gg <- vis.clonal.space(currHomeo_mat) +
##         ggtitle(paste0(currImmuno_v, " Clonal Space Homeostasis")) +
##         theme(plot.title = element_text(hjust = 0.5))

##     print(currHomeo_gg)
## } # for i

## rownames(immunoHomeo_mat) <- immunoRows_v

## immunoHomeo_gg <- vis.clonal.space(immunoHomeo_mat) +
##     ggtitle("Immunotherapy Clonal Space Homeostasis") +
##     theme(plot.title = element_text(hjust = 0.5))

## print(immunoHomeo_gg)

## dev.off()


## fullHomeo_mat <- clonal.space.homeostasis(.data = cloneData_lsdt, .prop.col = "Normalized clone fraction")

## nullSamples_v <- paste0("S", metadata_dt[immuno == "null", `sampleID/barcode`])
## aCD40Samples_v <- paste0("S", metadata_dt[immuno == "aCD40", `sampleID/barcode`])

## nullHomeo_mat <- fullHomeo_mat[rownames(fullHomeo_mat) %in% nullSamples_v,]
## aCD40Homeo_mat <- fullHomeo_mat[rownames(fullHomeo_mat) %in% aCD40Samples_v,]


## compareRare <- t.test(nullHomeo_mat[,1], aCD40Homeo_mat[,1], alternative = "two.sided")
## compareSmall <- t.test(nullHomeo_mat[,2], aCD40Homeo_mat[,2], alternative = "less")
## compareMedium <- t.test(nullHomeo_mat[,3], aCD40Homeo_mat[,3], alternative = "greater")
## compareLarge <- t.test(nullHomeo_mat[,4], aCD40Homeo_mat[,4], alternative = "greater")
## compareHyper <- t.test(nullHomeo_mat[,5], aCD40Homeo_mat[,5], alternative = "less")
