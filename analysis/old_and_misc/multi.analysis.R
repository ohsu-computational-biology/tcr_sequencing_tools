#	RScript that calculates number of unique clonotypes, Shannon entropy, and clonality for multiple clonotype files
#
#   Input is a directory containing one or more files with MiXCR output format
#   
#   Load necessary libraries
library(entropy);
library(data.table);

#	Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);
#   Directory should contain all and only .txt files from exportClones, post-normalization with the following
#   format:
#     Clone count   Clone fraction    Clonal sequence(s)    AA. Seq. CDR3   Best V Hit    Best J Hit    V segments
#     J segments    Normalized clone count    Normalized clone fraction

clone.dir <- arguments[1];    # Typically .../dhaarini/DNAXXXXLC/normalization/normalized_clones/
count.dir <- arguments[2];    # Typically .../dhaarini/DNAXXXXLC/spike_counts/9bp/counts/
out.dir <- arguments[3];
divisions_v <- arguments[4]; # How many top-clone divisions to make? Should be a comma-separated list of integers with no quotes nor spaces
                                        # example - 10,25,50,100




#	Examine the current directory for the files to process and sort them.
clone.files.in.dir <- list.files(clone.dir);
clone.files.in.dir <- clone.files.in.dir[order(as.numeric(gsub(".*_S|_.*$|\\..*$", '', clone.files.in.dir)))]
count.files.in.dir <- list.files(count.dir);
count.files.in.dir <- count.files.in.dir[order(as.numeric(gsub(".*_S|\\..*", '', count.files.in.dir)))]

if (length(clone.files.in.dir) != length(count.files.in.dir)) stop("Files do not match")


# Create empty arrays
raw.entropies <- numeric(length(clone.files.in.dir)); naive.entropies <- raw.entropies; nb.entropies <- raw.entropies
unique.clones <- NULL
raw.clonality <- NULL; naive.clonality <- NULL; nb.clonality <- NULL
raw.max.clonal.freq <- numeric(length(count.files.in.dir)); naive.max.clonal.freq <- raw.max.clonal.freq; nb.max.clonal.freq <- raw.max.clonal.freq
raw.norm.entropy <- numeric(length(clone.files.in.dir)); naive.norm.entropy <- raw.norm.entropy; nb.norm.entropy <- raw.norm.entropy
raw.adaptive.clonality <- NULL; naive.adaptive.clonality <- raw.adaptive.clonality; nb.adaptive.clonality <- raw.adaptive.clonality
raw.max.clone.count <- numeric(length(clone.files.in.dir)); naive.max.clone.count <- raw.max.clone.count; nb.max.clone.count <- raw.max.clone.count


###
### Handle different divisions
###
topDivisions_lsdf <- list()

### Split divisions into individual numbers
divisions_v <- sapply(strsplit(divisions_v, split = ',')[[1]], function(x) as.numeric(x), USE.NAMES=F)

### Create a list of list of data.frames, one list for each division containing 3 data.frames (1 for each norm type)
topDivisions_lslsdf <- as.list(divisions_v)
names(topDivisions_lslsdf) <- sapply(divisions_v, as.character)

lapply(topDivisions_lslsdf, function(x) {
    raw <- data.frame(matrix(nrow = x, ncol = length(clone.files.in.dir)))
    naive <- raw; nb <-	raw
    x <- list("raw" = raw, "naive" = naive, "nb" = nb)
    return(x)
})
## raw.top.10 <- data.frame(matrix(nrow = 10, ncol = length(clone.files.in.dir))); naive.top.10 <- raw.top.10; nb.top.10 <- raw.top.10
## raw.top.25 <- data.frame(matrix(nrow = 25, ncol = length(clone.files.in.dir))); naive.top.25 <- raw.top.25; nb.top.25 <- raw.top.25
## raw.top.50 <- data.frame(matrix(nrow = 50, ncol = length(clone.files.in.dir))); naive.top.50 <- raw.top.50; nb.top.50 <- raw.top.50
## raw.top.100 <- data.frame(matrix(nrow = 100, ncol = length(clone.files.in.dir))); naive.top.100 <- raw.top.100; nb.top.100 <- raw.top.100
## raw.top.200 <- data.frame(matrix(nrow = 200, ncol = length(clone.files.in.dir))); naive.top.200 <- raw.top.200; nb.top.200 <- raw.top.200
## raw.top.250 <- data.frame(matrix(nrow = 250, ncol = length(clone.files.in.dir))); naive.top.250 <- raw.top.250; nb.top.250 <- raw.top.250
## raw.top.500 <- data.frame(matrix(nrow = 500, ncol = length(clone.files.in.dir))); naive.top.500 <- raw.top.500; nb.top.500 <- raw.top.500

for(i in 1:length(clone.files.in.dir))	{

    ###
    ### DATA
    ###
    #   get a clone file to process
    clone.curr.file <- clone.files.in.dir[i];
    clone.index <- strsplit(clone.curr.file, split = "_|\\.")[[1]][2]

    clone.curr.record <- fread(file.path(clone.dir, clone.curr.file),
                            check.names=FALSE,
                            stringsAsFactors=FALSE);
    #   get a count file to process
    count.curr.file <- count.files.in.dir[i];
    count.index <- strsplit(count.curr.file, split = "_|\\.")[[1]][2]

    count.curr.record <- fread(file.path(count.dir, count.curr.file))

    if (count.index != clone.index) stop("Mismatching files, check index: ", i)

    # Depending on if original or data_subset, we need a norm fraction column
    if ("New.norm.fraction" %in% colnames(clone.curr.record)){
      column <- "New.norm.fraction"
    } else {
      column <- "Normalized clone fraction"
    } # if

 #    column <- "Normalized.clone.fraction"

    ###
    ### CALCULATIONS
    ###

    ## count number of lines in file, i.e. number of unique clonotypes
    unique.clones[i] <- clone.curr.record[,.N]

    ##   calculate entropy
    raw.entropies[i] <- entropy(clone.curr.record[["Clone fraction"]], method="ML", unit="log");
    naive.entropies[i] <- entropy(clone.curr.record[[column]], method="ML", unit="log");
    nb.entropies[i] <- entropy(clone.curr.record[["nb.clone.fraction"]], method="ML", unit="log");
    

    ##   calculate clonality
    raw.clonality[i] <- 1 - (raw.entropies[i] / log(unique.clones[i]))
    naive.clonality[i] <- 1 - (naive.entropies[i] / log(unique.clones[i]))
    nb.clonality[i] <- 1 - (nb.entropies[i] / log(unique.clones[i]))

    ##   calculate clonality as the inverse of normalized shannon entropy
    ## Normalized shannon entropy
    
    raw.norm.entropy[i] <- raw.entropies[i] / log(unique.clones[i])
    naive.norm.entropy[i] <- naive.entropies[i] / log(unique.clones[i])
    nb.norm.entropy[i] <- nb.entropies[i] / log(unique.clones[i])
    ## New clonality
    raw.adaptive.clonality[i] <- 1 / raw.norm.entropy[i]
    naive.adaptive.clonality[i] <- 1 / naive.norm.entropy[i]
    nb.adaptive.clonality[i] <- 1 / nb.norm.entropy[i]
	

    ## Change clone frequency column to a percentage
    clone.curr.record[["Clone fraction"]] <- clone.curr.record[["Clone fraction"]] * 100
    clone.curr.record[[column]] <- clone.curr.record[[column]] * 100
    clone.curr.record[["nb.clone.fraction"]] <- clone.curr.record[["nb.clone.fraction"]] * 100
    
    
    ##	Calculate Max. clonotype frequency
    raw.max.clonal.freq[i] <- round(max(clone.curr.record[["Clone fraction"]]), digits = 4)
    naive.max.clonal.freq[i] <- round(max(clone.curr.record[[column]]), digits = 4)
    nb.max.clonal.freq[i] <- round(max(clone.curr.record[["nb.clone.fraction"]]), digits = 4)
    

    ##   Record maximum clone count
    raw.max.clone.count[i] <- max(clone.curr.record$`Clone count`)
    naive.max.clone.count[i] <- max(clone.curr.record$`Normalized clone count`)
    nb.max.clone.count[i] <- max(clone.curr.record$`nb.clone.count`)

   #  Record frequencies for top 10 and top 25 clones
    clone.curr.record <- clone.curr.record[order(clone.curr.record[["Clone fraction"]], decreasing = T),]
    raw.top.10[,i] <- clone.curr.record[["Clone fraction"]][1:10]
    raw.top.25[,i] <- clone.curr.record[["Clone fraction"]][1:25]
    raw.top.50[,i] <- clone.curr.record[["Clone fraction"]][1:50]
    raw.top.100[,i] <- clone.curr.record[["Clone fraction"]][1:100]
    raw.top.200[,i] <- clone.curr.record[["Clone fraction"]][1:200]
    raw.top.250[,i] <- clone.curr.record[["Clone fraction"]][1:250]
    raw.top.500[,i] <- clone.curr.record[["Clone fraction"]][1:500]

    clone.curr.record <- clone.curr.record[order(clone.curr.record[[column]], decreasing = T),]
    naive.top.10[,i] <- clone.curr.record[[column]][1:10]
    naive.top.25[,i] <- clone.curr.record[[column]][1:25]
    naive.top.50[,i] <- clone.curr.record[[column]][1:50]
    naive.top.100[,i] <- clone.curr.record[[column]][1:100]
    naive.top.200[,i] <- clone.curr.record[[column]][1:200]
    naive.top.250[,i] <- clone.curr.record[[column]][1:250]
    naive.top.500[,i] <- clone.curr.record[[column]][1:500]

    clone.curr.record <- clone.curr.record[order(clone.curr.record[["nb.clone.fraction"]], decreasing = T),]
    nb.top.10[,i] <- clone.curr.record[["nb.clone.fraction"]][1:10]
    nb.top.25[,i] <- clone.curr.record[["nb.clone.fraction"]][1:25]
    nb.top.50[,i] <- clone.curr.record[["nb.clone.fraction"]][1:50]
    nb.top.100[,i] <- clone.curr.record[["nb.clone.fraction"]][1:100]
    nb.top.200[,i] <- clone.curr.record[["nb.clone.fraction"]][1:200]
    nb.top.250[,i] <- clone.curr.record[["nb.clone.fraction"]][1:250]
    nb.top.500[,i] <- clone.curr.record[["nb.clone.fraction"]][1:500]


    #   update progress
    if((i %%10) == 0)   {
        cat("Processing file ", i, " (out of ", length(count.files.in.dir), ")\n", sep="");
    }   #   fi

}	#	for i

# Summarize top10 and top25 data
raw.top.10.summary <- round(t(apply(raw.top.10, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
raw.top.25.summary <- round(t(apply(raw.top.25, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
raw.top.50.summary <- round(t(apply(raw.top.50, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
raw.top.100.summary <- round(t(apply(raw.top.100, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
raw.top.200.summary <- round(t(apply(raw.top.200, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
raw.top.250.summary <- round(t(apply(raw.top.250, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
raw.top.500.summary <- round(t(apply(raw.top.500, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)

naive.top.10.summary <- round(t(apply(naive.top.10, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
naive.top.25.summary <- round(t(apply(naive.top.25, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
naive.top.50.summary <- round(t(apply(naive.top.50, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
naive.top.100.summary <- round(t(apply(naive.top.100, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
naive.top.200.summary <- round(t(apply(naive.top.200, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
naive.top.250.summary <- round(t(apply(naive.top.250, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
naive.top.500.summary <- round(t(apply(naive.top.500, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)

nb.top.10.summary <- round(t(apply(nb.top.10, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
nb.top.25.summary <- round(t(apply(nb.top.25, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
nb.top.50.summary <- round(t(apply(nb.top.50, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
nb.top.100.summary <- round(t(apply(nb.top.100, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
nb.top.200.summary <- round(t(apply(nb.top.200, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
nb.top.250.summary <- round(t(apply(nb.top.250, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
nb.top.500.summary <- round(t(apply(nb.top.500, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)

#   create output data.frame
output.df <- data.frame(clone.files.in.dir, unique.clones, raw.entropies, naive.entropies, nb.entropies,
                        raw.norm.entropy, naive.norm.entropy, nb.norm.entropy,
                        raw.clonality, naive.clonality, nb.clonality,
                        raw.adaptive.clonality, naive.adaptive.clonality, nb.adaptive.clonality,
	     		raw.max.clonal.freq, naive.max.clonal.freq, nb.max.clonal.freq,
                        raw.max.clone.count, naive.max.clone.count, nb.max.clone.count,
			raw.top.10.summary, raw.top.25.summary, raw.top.50.summary, raw.top.100.summary,
                        raw.top.200.summary, raw.top.250.summary, raw.top.500.summary,
			naive.top.10.summary, naive.top.25.summary, naive.top.50.summary, naive.top.100.summary,
                        naive.top.200.summary, naive.top.250.summary, naive.top.500.summary,
			nb.top.10.summary, nb.top.25.summary, nb.top.50.summary, nb.top.100.summary,
                        nb.top.200.summary, nb.top.250.summary, nb.top.500.summary);

colnames(output.df) <- c("File", "Unique Clonotypes", "raw.shannon", "naive.shannon", "nb.shannon",
                         "raw.norm.shannon", "naive.norm.shannon", "nb.norm.shannon",
                         "raw.Clonality", "naive.Clonality", "nb.Clonality",
                         "raw.Adaptive.Clonality", "naive.Adaptive.Clonality", "nb.Adaptive.Clonality",
                         "raw.Max.Clonal.Freq", "naive.Max.Clonal.Freq", "nb.Max.Clonal.Freq",
                         "raw.Max.Clone.Count", "naive.Max.Clone.Count", "nb.Max.Clone.Count",
                         "raw.top10.mean", "raw.top10.median", "raw.top10.sum", "raw.top25.mean", "raw.top25.median", "raw.top25.sum",
                         "raw.top50.mean", "raw.top50.median", "raw.top50.sum", "raw.top100.mean", "raw.top100.median", "raw.top100.sum",
                         "raw.top200.mean", "raw.top200.median", "raw.top200.sum", "raw.top250.mean", "raw.top250.median", "raw.top250.sum",
                         "raw.top500.mean", "raw.top500.median", "raw.top500.sum",
                         "naive.top10.mean", "naive.top10.median", "naive.top10.sum", "naive.top25.mean", "naive.top25.median", "naive.top25.sum",
                         "naive.top50.mean", "naive.top50.median", "naive.top50.sum", "naive.top100.mean", "naive.top100.median", "naive.top100.sum",
                         "naive.top200.mean", "naive.top200.median", "naive.top200.sum", "naive.top250.mean", "naive.top250.median", "naive.top250.sum",
                         "naive.top500.mean", "naive.top500.median", "naive.top500.sum",
                         "nb.top10.mean", "nb.top10.median", "nb.top10.sum", "nb.top25.mean", "nb.top25.median", "nb.top25.sum",
                         "nb.top50.mean", "nb.top50.median", "nb.top50.sum", "nb.top100.mean", "nb.top100.median", "nb.top100.sum",
                         "nb.top200.mean", "nb.top200.median", "nb.top200.sum", "nb.top250.mean", "nb.top250.median", "nb.top250.sum",
                         "nb.top500.mean", "nb.top500.median", "nb.top500.sum")

#   write output
file.name <- "multi.uniques.shannon.clonality.txt"
write.table(output.df, 
            file=paste(out.dir, file.name, sep = ''),
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
