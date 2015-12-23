#	RScript that aggregates reports from MiXCR's alignment tool
#   
#   At this point the script simply accumulates results, but it'd be easy to add
#       some visualization, analysis, etc. once the data is aggregated

#   Load required libraries
library(stringr);

#	Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);
#   Directory should contain all and only MiXCR alignment reports for a given "batch"
working.dir <- arguments[1];

#   for debugging purposes
#working.dir <- "/Users/leyshock/Desktop/TCRseq/tools/temp/QC/";

#	Examine the current directory for the files to process
files.in.dir <- list.files(working.dir);

output.df <- data.frame(analysis.date=character(),  #   1
                        inputs=character(),
                        output=character(),
                        command=character(),
                        total.reads=integer(),  #   5
                        aligned.reads=integer(),
                        aligned.pct=numeric(),
                        filtered.diff.v.j.loci=numeric(),
                        failed.alignment.v.hits=numeric(),
                        failed.alignment.j.hits=numeric(),  #   10
                        failed.alignment.low.score=numeric(),
                        overlapped.pct=numeric(),
                        overlapped.and.aligned.pct=numeric(),
                        overlapped.and.not.aligned.pct=numeric(),   #   14
                        stringsAsFactors=FALSE);

for(i in 1:length(files.in.dir))	{
    #   get a QC file to process
    curr.file <- file.path(working.dir, files.in.dir[i]);
    
	curr.record <- readLines(curr.file);

    #   misc QC
    if(length(curr.record) != 15)   {
        stop("Unexpected length of report for file: ", curr.file, "\n", sep="");
    }   #   fi

    curr.date <- str_split(curr.record[1], ":")[[1]];
    curr.date <- curr.date[-1];
    curr.date <- paste(curr.date, collapse="");
    output.df[i,]$analysis.date <- curr.date;
    
    curr.input <- str_split(curr.record[2], ":")[[1]][2];
    curr.input <- str_trim(curr.input);
    curr.input <- basename(curr.input);
    output.df[i,]$inputs <- curr.input;

    curr.output <- str_split(curr.record[3], ":")[[1]][2];
    curr.output <- str_trim(curr.output);
    curr.output <- basename(curr.output);
    output.df[i,]$output <- curr.output;

    curr.command <- str_split(curr.record[4], ":")[[1]][2];
    curr.command <- str_trim(curr.command);
    curr.command <- basename(curr.command);
    output.df[i,]$command <- curr.command;
    
    curr.total <- str_split(curr.record[5], ":")[[1]][2];
    output.df[i,]$total.reads <- curr.total;
    
    curr.aligned <- str_split(curr.record[6], ":")[[1]][2];
    output.df[i,]$aligned.reads <- curr.aligned;

    for(j in 7:14)  {
        curr.temp <- str_split(curr.record[j], ":")[[1]][2];
        curr.temp <- as.numeric(str_replace(curr.temp, "%", ""));
        output.df[i,j]<- curr.temp;
    }   #   for j
     
}	#	for i

#   order results, to make them more intuitive
output.df <- output.df[order(output.df$aligned.pct, decreasing=TRUE),];

write.table(output.df, 
            file="mixcr.alignment.QC.summary.txt",
            quote=FALSE,
            sep=",",
            row.names=FALSE)

