#	RScript that aggregates reports from MiXCR's assembly tool
#   
#   At this point the script simply accumulates results, but it'd be easy to add
#       some visualization, analysis, etc. once the data is aggregated

#   Load required libraries
library(stringr);

#	Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);
#   Directory should contain all and only MiXCR alignment reports for a given "batch"
working.dir <- arguments[1]; # $data/mixcr/reports/assemble
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
                        clonotype.count=integer(),  
			avg.reads.per.clonotype=integer(),
                        num.reads.used=numeric(), pct.used.of.total=numeric(),
			num.reads.used.b4.clust=numeric(), pct.of.total=numeric(),	# 10
                        num.reads.used.as.core=numeric(), pct.of.used=numeric(),
                        num.reads.mapped.low.quality=numeric(), pct.mapped.of.used=numeric(),
                        num.PCR.error.clust=numeric(), pct.PCR.clust.of.used=numeric(),
			num.VJC.clust=numeric(), pct.VJC.clust.of.used=numeric(),
                        num.dropped.no.clonal.seq=numeric(), pct.dropped.no.clonal=numeric(),	# 15
                        num.reads.dropped.low.qual=numeric(), pct.dropped.low.quality=numeric(),
                        num.reads.dropped.fail.map=numeric(), pct.dropped.fail.map=numeric(),
			num.reads.drop.low.qual.clone=numeric(), pct.dropped.low.qual.clone=numeric(),
			clonotypes.elim.error.corr=numeric(),
			clonotypes.dropped.low.qual=numeric(),	# 20
			clonotypes.pre.clust.similar.VJC=numeric(),
                        stringsAsFactors=FALSE);

for(i in 1:length(files.in.dir))	{
    #   get a QC file to process
    curr.file <- file.path(working.dir, files.in.dir[i]);

    curr.record <- readLines(curr.file);

    #   misc QC
    if(length(curr.record) != 22)   {
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

    curr.command <- str_split(curr.record[6], ":")[[1]][2]
    curr.command <- str_trim(curr.command)
    output.df[i,]$command <- curr.command
    
    curr.count <- str_split(curr.record[7], ":")[[1]][2];
    output.df[i,]$clonotype.count <- curr.count;

    curr.avg.per.clone <- str_split(curr.record[8], ":")[[1]][2];
    output.df[i,]$avg.reads.per.clonotype <- curr.avg.per.clone;
    
#    curr.total <- str_split(curr.record[7], ":")[[1]][2];
#    output.df[i,]$total.reads.used <- curr.total;

    counter <- 0
    for(j in 9:18)  {
        curr.temp <- str_trim(str_split(curr.record[j], ":")[[1]][2]);
	curr.temp <- unlist(str_split(curr.temp, " "))
	curr.num <- as.numeric(curr.temp[1])
	curr.pct <- curr.temp[2]
        curr.pct <- as.numeric(str_replace_all(curr.pct, "\\(|%|\\)", ""));
        output.df[i,(j+counter)]<- curr.num;
	output.df[i,(j+counter+1)] <- curr.pct
	counter <- counter + 1
    }   #   for j

    for (j in 19:21) {
    	curr.temp <- str_split(curr.record[j], ":")[[1]][2];
	curr.temp <- as.numeric(str_trim(curr.temp));
	output.df[i,(j+counter)] <- curr.temp;
    } # for j
     
}	#	for i

#   order table to make results more intuitive
output.df <- output.df[order(output.df$pct.used.of.total, decreasing=TRUE),];

write.table(output.df, 
            file=file.path(out.dir, "mixcr.assemble.QC.summary.txt"),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

