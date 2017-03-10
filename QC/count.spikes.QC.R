#	RScript that aggregates QC outputs from count.spikes.R
#   
#   At this point the script simply accumulates results, but it'd be easy to add
#       some visualization, analysis, etc. once the data is aggregated

### Get command-line arguments
library(data.table)
arguments <- commandArgs(trailingOnly=TRUE);

working.dir <- arguments[1]; # $data/spike_counts/[bp]/qc/
out.dir <- arguments[2]


### Examine the current directory for the files to process
files.in.dir <- list.files(working.dir);

output.df <- data.frame();

for(i in 1:length(files.in.dir))	{
    ##   get a QC file to process
    curr.file <- paste(working.dir, files.in.dir[i], sep='');
    curr.record <- fread(curr.file)
    output.df <- rbind(output.df, curr.record);

}	#	for i

write.table(output.df, 
            file=file.path(out.dir, "count.spikes.QC.summary.txt"),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

