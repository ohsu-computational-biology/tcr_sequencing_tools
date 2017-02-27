#	RScript that aggregates QC outputs from count.spikes.R
#   
#   At this point the script simply accumulates results, but it'd be easy to add
#       some visualization, analysis, etc. once the data is aggregated

#	Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);
#   Directory should contain all and only QC files for a given "batch"
working.dir <- arguments[1];
out.dir <- arguments[2]

#   for debugging purposes
#working.dir <- "/Users/leyshock/Desktop/TCRseq/tools/temp/QC/";

#	Examine the current directory for the files to process
files.in.dir <- list.files(working.dir);

output.df <- data.frame();

for(i in 1:length(files.in.dir))	{
    #   get a QC file to process
    curr.file <- paste(working.dir, files.in.dir[i], sep='');

    curr.record <- read.table(curr.file,
                            header=TRUE,
                            sep=",",
                            stringsAsFactors=FALSE);
    output.df <- rbind(output.df, curr.record);

}	#	for i

write.table(output.df, 
            file=file.path(out.dir, "count.spikes.QC.summary.txt"),
            quote=FALSE,
            sep=",",
            row.names=FALSE)

