#   Script for normalizing R

#   disable scientific notation
options(scipen=999);

#   load required libraries
library(stringr);

qc.normalization <- function(path.to.raw.counts, path.to.normalized.counts) {

raw.counts <- list.files(path.to.raw.counts);
processed.counts <- list.files(path.to.normalized.counts);

#   check for parallelism of samples
sample.id.raw.counts <- character(length(raw.counts));
sample.id.processed.counts <- character(length(processed.counts));
for(i in 1:length(raw.counts))  {
    sample.id.raw.counts[i] <- str_split(raw.counts[i], "_")[[1]][1    ];
    sample.id.processed.counts[i] <- str_split(processed.counts[i], "_")[[1]][1];
}   #   for i

#   Check for parallelism of files
sample.comparison <- sample.id.raw.counts == sample.id.processed.counts;
sample.comparison <- which(sample.comparison == FALSE);
if(length(sample.comparison) > 0)   {
    stop("Mismatch between raw.counts and processed.counts\n");
}   #   fi

for(i in 1:length(raw.counts))  {
    #   Note that we use check.names=FALSE; this preserves the original column names,
    #       which is useful since some downstream tools (e.g. VDJTools' Convert() function)
    #       assume certain column names
    curr.raw <- read.table(file.path(path.to.raw.counts, raw.counts[i]),
                            check.names=FALSE,  
                            stringsAsFactors=FALSE,
                            sep="\t",
                            header=TRUE);

    curr.normalized <- read.table(file.path(path.to.normalized.counts, processed.counts[i]),
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
    combined.table$raw.file.name <- raw.counts[i];
    combined.table$normalized.file.name <- processed.counts[i];
    combined.table$normalization.factor <- round((combined.table$normalized.clone.count / combined.table$raw.clone.count), digits=1);

    output.file.name <- paste(sample.id.raw.counts[i], "_normalization_QC.txt", sep="");
    output.file.name <- file.path(path.to.normalized.counts, output.file.name);
    cat("Writing output to: ", output.file.name, "\n", sep="");

    write.table(combined.table,
                file=output.file.name,
                quote=FALSE,
                sep=",",
                row.names=FALSE);

    #   reset value
    rm(combined.table);
}   #   for i

}   #   qc.normalization()
