### Calculate Shannon Entropy and Clonality for top 50 Clones in a file

##   Input is a directory containing one or more files with MiXCR output format
   
##   Load necessary libraries
library(entropy);
library(data.table);

##	Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);

clone.dir <- arguments[1];    # Typically .../dhaarini/DNAXXXXLC/normalization/normalized_clones/
out.dir <- arguments[2];


##	Examine the current directory for the files to process and sort them.
clone.files.in.dir <- list.files(clone.dir);
clone.files.in.dir <- clone.files.in.dir[order(as.numeric(gsub(".*_S|_.*$|\\..*$", '', clone.files.in.dir)))]

## Create empty arrays
clone.indexes <- NULL
calculated.entropies <- numeric(length(clone.files.in.dir));
clonality <- numeric(length(clone.files.in.dir));
max.clonal.freq <- numeric(length(count.files.in.dir));
max.clone.count <- numeric(length(clone.files.in.dir));

for(i in 1:length(clone.files.in.dir))	{

    ##   get a clone file to process
    clone.curr.file <- clone.files.in.dir[i];
    clone.index <- strsplit(clone.curr.file, split = "_|\\.")[[1]][2]
    clone.indexes[i] <- clone.index
    clone.curr.record <- fread(file.path(clone.dir, clone.curr.file))
    
    ## Depending on if original or data_subset, we need a norm fraction column
    if ("New.norm.fraction" %in% colnames(clone.curr.record)){
        column <- "New.norm.fraction"
    } else {
        column <- "Normalized clone fraction"
    } # if

    ## Sort By Frequnecy and Subset for Top 50
    clone.curr.record <- clone.curr.record[order(clone.curr.record[[column]], decreasing = T),]
    top50.clones <- clone.curr.record[1:50,]

    ## Calculate Entropy
    calculated.entropies[i] <- entropy(top50.clones[[column]], method = "ML", unit = "log")

    ## Calculate Clonality
    clonality[i] <- 1 - (calculated.entropies[i] / log(50))
        
    ##	Calculate Max. clonotype frequency
    max.clonal.freq[i] <- round(max(top50.clones[[column]]), digits = 4)

    ##   Record maximum clone count
    max.clone.count[i] <- max(top50.clones$`Normalized clone count`)

    ##   update progress
    if((i %%10) == 0)   {
        cat("Processing file ", i, " (out of ", length(count.files.in.dir), ")\n", sep="");
    }   #   fi

}	#	for i

###   create output data.frame
output.df <- data.frame(clone.indexes, calculated.entropies, clonality, max.clonal.freq, max.clone.count)
colnames(output.df) <- c("File", "Shannon Entropy", "Clonality", "Max Clonal Freq", "Max Clone Count")

### Write Output
file.name <- "top50.analysis.txt"
write.table(output.df, 
            file=paste(out.dir, file.name, sep = ''),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
