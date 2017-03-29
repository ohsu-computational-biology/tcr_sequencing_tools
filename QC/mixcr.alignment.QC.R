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
out.dir <- arguments[2]

#   for debugging purposes
#working.dir <- "/Users/leyshock/Desktop/TCRseq/tools/temp/QC/";

#	Examine the current directory for the files to process
files.in.dir <- list.files(working.dir);

output.df <- data.frame(analysis.date=character(),  #   1
                        inputs=character(),
                        output=character(),
			version=character(),
			time=character(),	# 5
                        command=character(),
                        total.reads=integer(),  
                        aligned.reads=integer(), aligned.pct=numeric(),
                        failed.alignment.no.hits=numeric(), pct.no.hits=numeric(),
                        failed.alignment.no.j.hits=numeric(), pct.no.j.hits=numeric(),	# 10
                        failed.alignment.low.score=numeric(), pct.low.score=numeric(),
                        num.overlapped=numeric(), pct.overlapped=numeric(),
                        num.overlapped.and.aligned=numeric(), pct.overlapped.and.aligned=numeric(),
                        num.overlapped.and.not.aligned=numeric(), pct.overlapped.and.not.aligned=numeric(),   #   14
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

    curr.version <- str_split(curr.record[4], ":|;")[[1]][2];
    curr.version <- str_trim(curr.version);
    output.df[i,]$version <- curr.version;

    curr.time <- str_split(curr.record[5], ":")[[1]][2]
    curr.time <- str_trim(curr.time)
    output.df[i,]$time <- curr.time

    curr.command <- str_split(curr.record[6], ":")[[1]][2];
    curr.command <- str_trim(curr.command);
    ##    curr.command <- basename(curr.command);
    curr.command <- gsub(" --report.*", "", curr.command)
    output.df[i,]$command <- curr.command;
    
    curr.total <- str_split(curr.record[7], ":")[[1]][2];
    output.df[i,]$total.reads <- curr.total;
    
#    curr.aligned <- str_split(curr.record[6], ":")[[1]][2];
#    output.df[i,]$aligned.reads <- curr.aligned;

    counter <- 0
    for(j in 8:14)  {
        curr.temp <- str_trim(str_split(curr.record[j], ":")[[1]][2]);
	curr.temp <- unlist(str_split(curr.temp, " "));
	curr.num <- as.numeric(curr.temp[1])
	curr.pct <- curr.temp[2]
	curr.pct <- as.numeric(str_replace_all(curr.pct, "\\(|%|\\)", ""))
        output.df[i,(j+counter)]<- curr.num;
	output.df[i,(j+counter+1)] <- curr.pct
	counter <- counter + 1
    }   #   for j
     
}	#	for i

#   order results, to make them more intuitive
output.df <- output.df[order(output.df$aligned.pct, decreasing=TRUE),];

write.table(output.df, 
            file=file.path(out.dir, "mixcr.alignment.QC.summary.txt"),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

