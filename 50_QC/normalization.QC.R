#!/usr/bin/Rscript
###   Rscript that performs quality control on the normalization step of the TCRseq workflow

####################
### DEPENDENCIES ###
####################

suppressMessages(library(optparse))
suppressMessages(library(data.table))
options(scipen=999); # disable scientific notation

####################
### COMMAND ARGS ###
####################

optlist <- list(
        make_option(
                c("-r", "--rawDir"),
                type = "character",
                help = "Path to directory containing all and only raw clone count files for a batch."
        ),
	make_option(
		c("-n", "--normDir"),
		type = "character",
		help = "Path to directory containing all and only normalized clone count files for a batch."
	),
        make_option(
                c("-o", "--outDir"),
                type = "character",
                help = "Path to directory where results will be written."
        )
)

##################
### PARSE ARGS ###
##################

p <- OptionParser(usage = "%prog -r rawDir -n normDir -o outDir",
                option_list = optlist)
args <- parse_args(p)
opt <- args$options

rawDir_v <- args$rawDir
normDir_v <- args$normDir
outDir_v <- args$outDir

###############
### ACTIONS ###
###############
#arguments <- commandArgs(trailingOnly=TRUE);
#path.to.raw.clone.counts <- arguments[1];
#path.to.normalized.clone.counts <- arguments[2];
#out.dir <- arguments[3]

### Get files
rawFiles_v <- list.files(rawDir_v)
normFiles_v <- list.files(normDir_v)

### Sort
rawFiles_v <- rawFiles_v[order(as.numeric(gsub("^.*_S|_[a-z]+|\\.clono.*|.txt", "", rawFiles_v)))]
normFiles_v <- normFiles_v[order(as.numeric(gsub("^.*_S|_[a-z]+|\\.clono.*|.txt", "", normFiles_v)))]

### Get IDs
rawID_v <- gsub("^.*_S|_[a-z]+|.txt", "", rawFiles_v)
normID_v <- gsub("^.*_S|_[a-z]+|.txt", "", normFiles_v)

### Check parallelism
notMatched_v <- which(rawID_v != normID_v)
if (length(notMatched_v) > 0) stop("Mismatch between raw counts and norm counts")

### Create final output
output_df <- data.frame()

### Iterate over each sample
for(i in 1:length(rawFiles_v))  {

    ## Get files
    currRaw_v <- rawFiles_v[i]
    currNorm_v <- normFiles_v[i]

    ## Read files
    currRaw_dt <- fread(file.path(rawDir_v, currRaw_v))
    currNorm_dt <- fread(file.path(normDir_v, currNorm_v))

    ## Get column names
    cdr3Col_v <- grep("AA. Seq. CDR3|aaSeqCDR3", colnames(currRaw_dt), value = T)
    rawCountCol_v <- grep("Clone count|cloneCount", colnames(currRaw_dt), value = T)
    rawFreqCol_v <- grep("Clone fraction|cloneFraction", colnames(currRaw_dt), value = T)
    normCountCol_v <- grep("Normalized clone count", colnames(currNorm_dt), value = T)
    normFreqCol_v <- grep("Normalized clone fraction", colnames(currNorm_dt), value = T)
    nbCountCol_v <- grep("nb.clone.count", colnames(currNorm_dt), value = T)
    nbFreqCol_v <- grep("nb.clone.fraction", colnames(currNorm_dt), value = T)

    ## Basic QC - check that CDR3 sequences are teh same
    currRaw_CDR3 <- currRaw_dt[[cdr3Col_v]]
    currNorm_CDR3 <- currNorm_dt[[cdr3Col_v]]
    
    if(!identical(currRaw_CDR3, currNorm_CDR3)) {
        stop("Mistmatch between amino acid CDR3 region, raw and normalized");
    }   #   fi

    ## Get raw count and freq
    combinedTable_df <- data.frame(raw.clone.count = currRaw_dt[[rawCountCol_v]])
    combinedTable_df$raw.clone.percent <- currRaw_dt[[rawFreqCol_v]]

    ## Create variable to hold norm factor
    normFactorSummary_v <- NULL

    ## Get original norm method count and freq (if exists) also norm factor
    if (length(normCountCol_v) > 0) {
	combinedTable_df$medianNorm.count <- currNorm_dt[[normCountCol_v]]
	combinedTable_df$medianNorm.freq <- currNorm_dt[[normFreqCol_v]]
	medianNormFactor_v <- round((combinedTable_df$medianNorm.count / 
							combinedTable_df$raw.clone.count), digits = 1)
	medianFactorSummary_v <- as.vector(summary(medianNormFactor_v))[c(1,3,6)]
	names(medianFactorSummary_v) <- c("medianNorm.factor.min", "medianNorm.factor.median", "medianNorm.factor.max")
	normFactorSummary_v <- c(normFactorSummary_v, medianFactorSummary_v)
    }

    ## Get nb norm method count and freq (if exists) also norm factor
    if (length(nbCountCol_v) > 0){
        combinedTable_df$nbNorm.count <- currNorm_dt[[nbCountCol_v]]
        combinedTable_df$nbNorm.freq <- currNorm_dt[[nbFreqCol_v]]
        nbNormFactor_v <- round((combinedTable_df$nbNorm.count /
                                                        combinedTable_df$raw.clone.count), digits = 1)
	nbFactorSummary_v <- as.vector(summary(nbNormFactor_v))[c(1,3,6)]
	names(nbFactorSummary_v) <- c("nbNorm.factor.min", "nbNorm.factor.median", "nbNorm.factor.max")
	normFactorSummary_v <- c(normFactorSummary_v, nbFactorSummary_v)
    }

    ## Re-order
    if ((length(normCountCol_v) > 0) & (length(nbCountCol_v) > 0)) {
	combinedTable_df <- combinedTable_df[,c(1,3,5,2,4,6)]
    } else if ( ((length(normCountCol_v) > 0) & (length(nbCountCol_v) == 0))  |
		((length(normCountCol_v) == 0) & (length(nbCountCol_v) > 0)) ) {
	combinedTable_df <- combinedTable_df[,c(1,3,2,4)]
    } else if ((length(normCountCol_v) == 0) & (length(nbCountCol_v) == 0)) {
	stop("Neither 'Median Norm' nor 'NB Norm' columns are set. Please run normalization before using this QC tool.")
    } # fi

    ## Summarize
    combinedTable_v <- apply(combinedTable_df, 2, mean)
    combinedTable_v <- c(combinedTable_v, normFactorSummary_v)

    ## Extract column names
    colNames_v <- names(combinedTable_v)

    ## Add to final output
    output_df <- rbind(output_df, combinedTable_v)
    colnames(output_df) <- colNames_v

}   #   for i

### Add sample numbers
output_df <- cbind(rawID_v, output_df)

write.table(output_df,
	file = file.path(outDir_v, "normFactorQC.txt"),
	quote = F, row.names = F, sep = '\t')
