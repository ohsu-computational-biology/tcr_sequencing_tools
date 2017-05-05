#   Rscript that performs quality control on the normalization step of the
#       TCRseq workflow
#
#   Run this as an Rscript, with two arguments
#
#   Output is a csv file, suitable for further processing, or import into Excel

#   disable scientific notation
options(scipen=999);

###   load required libraries
library(data.table)
#library(stringr);

arguments <- commandArgs(trailingOnly=TRUE);
path.to.raw.clone.counts <- arguments[1];
path.to.normalized.clone.counts <- arguments[2];
out.dir <- arguments[3]


raw.clone.counts <- list.files(path.to.raw.clone.counts);
processed.clone.counts <- list.files(path.to.normalized.clone.counts);

###   check for parallelism of samples
sample.id.raw.clone.counts <- character(length(raw.clone.counts));
sample.id.processed.clone.counts <- character(length(processed.clone.counts));

###	TODO:  fix this, it relies on file naming convention
sample.id.raw.clone.counts <- sapply(raw.clone.counts, function(x) strsplit(x, "_")[[1]][2], USE.NAMES = F)
sample.id.processed.clone.counts <- sapply(processed.clone.counts, function(x) strsplit(x, "_")[[1]][2], USE.NAMES = F)

###   Check for parallelism of files
sample.comparison <- sample.id.raw.clone.counts == sample.id.processed.clone.counts;
sample.comparison <- which(sample.comparison == FALSE);
if(length(sample.comparison) > 0)   {
    stop("Mismatch between raw.clone.counts and processed.clone.counts\n");
}   #   fi

for(i in 1:length(raw.clone.counts))  {
    ## Read files
    curr.raw <- fread(file.path(path.to.raw.clone.counts, raw.clone.counts[i]))
    curr.normalized <- fread(file.path(path.to.normalized.clone.counts, processed.clone.counts[i]))

    ##   Basic QC
    curr.raw.CDR3 <- curr.raw$"AA. seq. CDR3";
    curr.normalized.CDR3 <- curr.normalized$"AA. seq. CDR3";
    
    if(!identical(curr.raw.CDR3, curr.normalized.CDR3)) {
        stop("Mistmatch between amino acid CDR3 region, raw and normalized");
    }   #   fi

    combined.table <- data.frame(raw.clone.count=curr.raw$"Clone count");
    combined.table$normalized.clone.count <- curr.normalized$"Normalized clone count";
    combined.table$raw.clone.percent <- curr.raw$"Clone fraction";
    combined.table$normalized.clone.percent <- curr.normalized$"Normalized clone fraction";
    combined.table$raw.file.name <- raw.clone.counts[i];
    combined.table$normalized.file.name <- processed.clone.counts[i];
    combined.table$normalization.factor <- round((combined.table$normalized.clone.count / combined.table$raw.clone.count), digits=1);

    output.file.name <- paste(sample.id.raw.clone.counts[i], "_normalization_QC.txt", sep="");
#    output.file.name <- file.path(path.to.normalized.clone.counts, output.file.name);
    cat("Writing output to: ", file.path(out.dir,output.file.name), "\n", sep="");

    write.table(combined.table,
                file=file.path(out.dir, output.file.name),
                quote=FALSE,
                sep="\t",
                row.names=FALSE);

    #   reset value
    rm(combined.table);
}   #   for i

