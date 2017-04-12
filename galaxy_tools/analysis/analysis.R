#	RScript that calculates number of unique clonotypes, Shannon entropy, and clonality for multiple clonotype files
#
#   Input is a directory containing one or more files with MiXCR output format
#   
#   Load necessary libraries
#source("https://cran.r-project.org/web/packages/entropy/index.html", echo=F, verbose=F)
suppressMessages(install.packages("entropy", repos = "http://cran.us.r-project.org", quiet = TRUE))
suppressMessages(library(entropy))
suppressMessages(library(data.table))

#	Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);
#   Directory should contain all and only .txt files from exportClones, post-normalization with the following
#   format:
#     Clone count   Clone fraction    Clonal sequence(s)    AA. Seq. CDR3   Best V Hit    Best J Hit    V segments
#     J segments    Normalized clone count    Normalized clone fraction

clone.dir <- arguments[1];    # Typically .../dhaarini/DNAXXXXLC/normalization/normalized_clones/
print(clone.dir)
clone.files.in.dir <- unlist(strsplit(clone.dir, ','))
print(clone.files.in.dir)

count.dir <- arguments[2];    # Typically .../dhaarini/DNAXXXXLC/spike_counts/9bp/counts/
print(count.dir)
count.files.in.dir <- unlist(strsplit(count.dir, ','))
print(count.files.in.dir)

output <- arguments[3];


# Create empty arrays
calculated.entropies <- NULL
unique.clones <- NULL
clonality <- NULL
max.clonal.freq <- NULL
norm.entropy <- NULL
adaptive.clonality <- NULL
max.clone.count <- NULL
top.10 <- data.frame(matrix(nrow = 10, ncol = length(clone.files.in.dir)))
top.25 <- data.frame(matrix(nrow = 25, ncol = length(clone.files.in.dir)))
top.50 <- data.frame(matrix(nrow = 50, ncol = length(clone.files.in.dir)))

for(i in 1:length(clone.files.in.dir))	{
    ##   get a clone file to process
    clone.curr.file <- clone.files.in.dir[i];

    clone.curr.record <- fread(clone.curr.file)
    
    ##   get a count file to process
    count.curr.file <- count.files.in.dir[i];

    count.curr.record <- fread(count.curr.file)

    ## Depending on if original or data_subset, we need a norm fraction column
    if ("New.norm.fraction" %in% colnames(clone.curr.record)){
        column <- "New.norm.fraction"
    } else {
        column <- "Normalized clone fraction"
    } # if

    ##
    ## Calculations
    ##
    
    unique.clones[i] <- clone.curr.record[,.N]

    ##   calculate entropy
    calculated.entropies[i] <- entropy(clone.curr.record[[column]], method="ML", unit="log");

    ##   calculate clonality
    clonality[i] <- 1 - (calculated.entropies[i] / log(unique.clones[i]))

    ## Calculate clonality as inverse of "normalized" shannon entropy
    ## Normalized shannon entropy
    norm.entropy[i] <- calculated.entropies[i] / log(unique.clones[i])
    ## Adaptive clonality
    adaptive.clonality[i] <- 1 / norm.entropy[i]

    ## Change clone freq column to a percentage
    clone.curr.record[[column]] <- clone.curr.record[[column]] * 100

    ##	Calculate Max. clonotype frequency
    max.clonal.freq[i] <- round(max(clone.curr.record[[column]]), digits = 4)

    ## Mac lone count
    max.clone.count[i] <- max(clone.curr.record$`Normalized clone count`)

    ##  Record frequencies for top 10 and top 25 clones
    clone.curr.record <- clone.curr.record[order(clone.curr.record[[column]], decreasing = T),]
    top.10[,i] <- clone.curr.record[[column]][1:10]
    top.25[,i] <- clone.curr.record[[column]][1:25]
    top.50[,i] <- clone.curr.record[[column]][1:50]

} # for i

### Summarize top10 and top25 data
top.10.summary <- round(t(apply(top.10, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
top.25.summary <- round(t(apply(top.25, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)
top.50.summary <- round(t(apply(top.50, 2, function(x) c(mean(x), median(x), sum(x)))), digits = 4)

## For debug
## cat("Clone Files\n", length(clone.files.in.dir)); cat("\n")
## cat("Entropies\n", length(calculated.entropies)); cat("\n")
## cat("Norm.Entropy\n", length(norm.entropy)); cat("\n")
## cat("unique.clones\n", length(unique.clones)); cat("\n")
## cat("Clonality\n", length(clonality)); cat("\n")
## cat("Adpative.clonality\n", length(adaptive.clonality)); cat("\n")
## cat("max.clonal.freq\n", length(max.clonal.freq)); cat("\n")
## cat("max.clone.count\n", length(max.clone.count)); cat("\n")
## cat("top.10.summary\n"); cat("\n")
## cat(length(top.10.summary)); cat("\n")
## cat("top.25.summary"); cat("\n")
## cat(length(top.25.summary)); cat("\n")
## cat("top.50.summary\n")
## cat(length(top.50.summary)); cat("\n")

###   create output data.frame
output.df <- data.frame(clone.files.in.dir, calculated.entropies, norm.entropy, unique.clones, clonality,
                        adaptive.clonality, max.clonal.freq, max.clone.count,
                        top.10.summary, top.25.summary, top.50.summary);

colnames(output.df) <- c("File", "Shannon Entropy", "Normalized Entropy", "Unique Clonotypes", "Clonality",
                         "Adaptive Clonality", "Max Clonal Freq", "Max Clone Count", "top10.mean", "top10.median",
                                                  "top10.sum", "top25.mean", "top25.median", "top25.sum", "top50.mean", "top50.median", "top50.sum");

###   write output
write.table(output.df, 
            file=output,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
