# Script to compare two biological replicates

# For DNA151124LC replicates are separated by 85 samples. S1 and S86, S2 and S87, etc.

# Get arguments
arguments <- commandArgs(trailingOnly = TRUE);

# Assign arguments
clones.dir <- arguments[1]    # Directory of normalized clone files (_alignment_clones_exported_normalized.txt)
counts.dir <- arguments[2]    # Directory of 9bp counts files (assembled.spike.counts.9bp.txt)
offset <- arguments[3]        # Number indicating distance between replicates
                              # ex) if 1-85 are replicate A and 86-170 are replicate B, offset = 85

# Create arrays of files from directories and sort them
clones.files <- list.files(clones.dir)
clones.files <- clones.files[order(as.numeric(gsub(".*_S|_alignment_.*", '', clones.files)))]
counts.files <- list.files(counts.dir)
counts.files <- counts.files[order(as.numeric(gsub(".*_S|\\..*", '', counts.files)))]

# Create empty arrays of appropriate length for calculations
sample1.uniques <- numeric(length(clones.files)/2);
sample2.uniques <- numeric(length(clones.files)/2);
tcr.counts.1 <- numeric(length(clones.files)/2);
tcr.counts.2 <- numeric(length(clones.files)/2);
tcr.counts.diff <- numeric(length(clones.files)/2)
difference <- numeric(length(clones.files)/2);
id <- numeric(length(clones.files)/2);
spike.counts.1 <- numeric(length(counts.files)/2);
spike.counts.2 <- numeric(length(counts.files)/2);
spike.diff <- numeric(length(counts.files)/2);


# Calculate various statistics for each pair of replicates
for (i in 1:(length(clones.files)/2)){
  curr.clone <- clones.files[i]
  rep.clone <- clones.files[i+offset]
  curr.count <- counts.files[i]
  rep.count <- counts.files[i+offset]
  curr.clone.data <- read.delim(file.path(clones.dir, curr.clone),
                                 check.names=FALSE,
                                 stringsAsFactors=FALSE);
  rep.clone.data <- read.delim(file.path(clones.dir, rep.clone),
                                check.names=FALSE,
                                stringsAsFactors=FALSE);
  curr.count.data <- read.delim(file.path(counts.dir, curr.count), sep = ',',
                          check.names=FALSE,
                          stringsAsFactors=FALSE);
  rep.count.data <- read.delim(file.path(counts.dir, rep.count), sep = ',',
                              check.names=FALSE,
                              stringsAsFactors=FALSE)
  curr.system.call.s1 <- paste("wc -l ", clones.dir, clones.files[i],
                               " | awk '{print $1}'", sep='');
  curr.system.call.s2 <- paste("wc -l ", clones.dir, clones.files[i+offset],
                               " | awk '{print $1}'", sep='');
  
  # Compare Unique Clonotypes
  sample1.uniques[i] <- as.numeric(system(curr.system.call.s1, intern=TRUE)) -1;
  sample2.uniques[i] <- as.numeric(system(curr.system.call.s2, intern=TRUE)) -1;
  difference[i] <- (sample1.uniques[i] - sample2.uniques[i]) / sample1.uniques[i] * 100
  
  # Compare Spike Counts
  spike.counts.1[i] <- curr.count.data[1,5]
  spike.counts.2[i] <- rep.count.data[1,5]
  spike.diff[i] <- (spike.counts.1[i] - spike.counts.2[i]) / spike.counts.1[i] * 100
  
  # Compare TCR Counts
  tcr.counts.1[i] <- sum(curr.clone.data$`Clone count`)
  tcr.counts.2[i] <- sum(rep.clone.data$`Clone count`)
  tcr.counts.diff[i] <- (tcr.counts.1[i] - tcr.counts.2[i]) / tcr.counts.1[i] * 100
  
  # Create ID array. Replicate1-Replicate2
  id[i] <- paste("S", i, "-", "S", i+offset, sep="");
}

# Create output data frame  
output.df <- data.frame(id, sample1.uniques, sample2.uniques, difference, spike.counts.1, spike.counts.2,
                        spike.diff, tcr.counts.1, tcr.counts.2, tcr.counts.diff);

#   write output
write.table(output.df, 
            file="replicate.comparison.txt",
            quote=FALSE,
            sep=",",
            row.names=FALSE)
