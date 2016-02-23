#   Rscript that performs quality control on the normalization step of the
#       TCRseq workflow
#
#   Run this as an Rscript, with two arguments
#
#   Output is a csv file, suitable for further processing, or import into Excel

#   disable scientific notation
options(scipen=999);

#   load required libraries
library(stringr);

arguments <- commandArgs(trailingOnly=TRUE);
path.to.raw.clone.counts <- arguments[1];
path.to.normalized.clone.counts <- arguments[2];


raw.clone.counts <- list.files(path.to.raw.clone.counts);
processed.clone.counts <- list.files(path.to.normalized.clone.counts);

#   check for parallelism of samples
sample.id.raw.clone.counts <- character(length(raw.clone.counts));
sample.id.processed.clone.counts <- character(length(processed.clone.counts));
#	TODO:  fix this, it relies on file naming convention
for(i in 1:length(raw.clone.counts))  {
    sample.id.raw.clone.counts[i] <- str_split(raw.clone.counts[i], "_")[[1]][2];
    sample.id.processed.clone.counts[i] <- str_split(processed.clone.counts[i], "_")[[1]][2];
}   #   for i

#   Check for parallelism of files
sample.comparison <- sample.id.raw.clone.counts == sample.id.processed.clone.counts;
sample.comparison <- which(sample.comparison == FALSE);
if(length(sample.comparison) > 0)   {
    stop("Mismatch between raw.clone.counts and processed.clone.counts\n");
}   #   fi

for(i in 1:length(raw.clone.counts))  {
    #   Note that we use check.names=FALSE; this preserves the original column names,
    #       which is useful since some downstream tools (e.g. VDJTools' Convert() function)
    #       assume certain column names
    curr.raw <- read.table(file.path(path.to.raw.clone.counts, raw.clone.counts[i]),
                            check.names=FALSE,  
                            stringsAsFactors=FALSE,
                            sep="\t",
                            header=TRUE);

    curr.normalized <- read.table(file.path(path.to.normalized.clone.counts, processed.clone.counts[i]),
                            check.names=FALSE,  
                            stringsAsFactors=FALSE,
                            sep="\t",
                            header=TRUE);

    #   Basic QC
    curr.raw.CDR3 <- curr.raw$"AA. seq. CDR3";
    curr.normalized.CDR3 <- curr.normalized$"AA. seq. CDR3";
    if(!identical(curr.raw.CDR3, curr.normalized.CDR3)) {
        stop("Mistmatch between amino acid CDR3 region, raw and normalized");
    }   #   fi

    combined.table <- data.frame(raw.clone.count=curr.raw$"Clone count");
    combined.table$normalized.clone.count <- curr.normalized$"Clone count";
    combined.table$raw.clone.percent <- curr.raw$"Clone fraction";
    combined.table$normalized.clone.percent <- curr.normalized$"Clone fraction";
    combined.table$raw.file.name <- raw.clone.counts[i];
    combined.table$normalized.file.name <- processed.clone.counts[i];
    combined.table$normalization.factor <- round((combined.table$normalized.clone.count / combined.table$raw.clone.count), digits=1);

    output.file.name <- paste(sample.id.raw.clone.counts[i], "_normalization_QC.txt", sep="");
    output.file.name <- file.path(path.to.normalized.clone.counts, output.file.name);
    cat("Writing output to: ", output.file.name, "\n", sep="");

    write.table(combined.table,
                file=output.file.name,
                quote=FALSE,
                sep=",",
                row.names=FALSE);

    #   reset value
    rm(combined.table);
}   #   for i

