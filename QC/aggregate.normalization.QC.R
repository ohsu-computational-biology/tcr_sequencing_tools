#   Script for normalizing R

#   disable scientific notation
options(scipen=999);

#   load required libraries
library(stringr);

qc.normalization <- function(path.to.QC.files)  {

all.files <- list.files(path.to.QC.files);

#   check for parallelism of samples
sample.ids <- character(length(all.files));
for(i in 1:length(sample.ids))  {
    sample.ids[i] <- str_split(all.files[i], "_")[[1]][1];
}   #   for i

output.data <- NULL;

for(i in 1:length(sample.ids))  {
    #   Note that we use check.names=FALSE; this preserves the original column names,
    #       which is useful since some downstream tools (e.g. VDJTools' Convert() function)
    #       assume certain column names
    curr.data <- read.table(file.path(path.to.QC.files, all.files[i]),
                            stringsAsFactors=FALSE,
                            sep=",",
                            header=TRUE);

    curr.normalization.factor <- curr.data$normalization.factor;
    output.data <- rbind(output.data, summary(curr.normalization.factor));

}   #   for i

output.data <- as.data.frame(output.data);
output.data$sample.id <- sample.ids;
#   arrange data so it's more intuitively presented
output.data <- output.data[order(output.data$Mean, decreasing=TRUE),];

output.file.name <- file.path(path.to.QC.files, "aggregate_normalization_factor_QC.txt");
cat("Writing output to: ", output.file.name, "\n", sep="");

write.table(output.data,
            file=output.file.name,
            quote=FALSE,
            sep=",",
            row.names=FALSE);

}   #   qc.normalization()
