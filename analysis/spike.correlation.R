# Script to calculate correlation of spike counts between samples

# Attempt at iteration
count.directory <- "~/Desktop/OHSU/TCR/tcr_spike/data/correlation/counts/"
offset <- 85

# Extract list of files from directory
count.files <- list.files(count.directory);

# Order them numerically
count.files <- count.files[order(as.numeric(gsub(".*_S|\\.assembled.*", '', count.files)))]


# Create empty data frame to populate with correlation coefficients
rep1 <- data.frame(matrix(ncol=85, nrow=85))
colnames(rep1) <- c(1:85)
rownames(rep1) <- c(1:85)
rep2 <- data.frame(matrix(ncol=85, nrow=85))
colnames(rep2) <- c(86:170)
rownames(rep2) <- c(86:170)

# Create matrix for first set of replicates
for (i in 1:(length(count.files)/2)){
  curr.count <- count.files[i]
  curr.count.data <- read.delim(file.path(count.directory, curr.count), sep = ',',
                                check.names=FALSE,
                                stringsAsFactors=FALSE);
  curr.spikes <- curr.count.data$spike.count
  for (j in 1:(length(count.files)/2)){
      curr.compare.count <- count.files[j]
      curr.compare.count.data <- read.delim(file.path(count.directory, curr.compare.count), sep = ',',
                                       check.names=FALSE, stringsAsFactors=FALSE);
      curr.compare.spikes <- curr.compare.count.data$spike.count
      rep1[i,j] <- summary(lm(curr.spikes ~ curr.compare.spikes))$r.squared
  }
}

# Create matrix for second set of replicates
for (i in 1:(length(count.files)/2)){
  curr.count <- count.files[i+offset]
  curr.count.data <- read.delim(file.path(count.directory, curr.count), sep = ',',
                                check.names=FALSE,
                                stringsAsFactors=FALSE);
  curr.spikes <- curr.count.data$spike.count
  for (j in 1:(length(count.files)/2)){
    #if (j < 85){
      curr.compare.count <- count.files[j+offset]
      curr.compare.count.data <- read.delim(file.path(count.directory, curr.compare.count), sep = ',',
                                        check.names=FALSE, stringsAsFactors=FALSE);
      curr.compare.spikes <- curr.compare.count.data$spike.count
      rep2[i,j] <- summary(lm(curr.spikes ~ curr.compare.spikes))$r.squared
    #} # fi
  } # for j
} # for i

# Calculate average correlation coefficient
# Samples 1-85
avg.r2.rep1 <- length(1:85)
for (i in 1:85){
  avg.r2.rep1[i] <- sum(rep2[i,])
}
avg.r2.rep1 <- sum(avg.r2)
avg.r2.rep1 <- avg.r2 / 7225

# Samples 86-170
avg.r2.rep2 <- length(1:85)
for (i in 1:85){
  avg.r2.rep2[i] <- sum(rep2[i,])
}
avg.r2.rep2 <- sum(avg.r2)
avg.r2.rep2 <- avg.r2 / 7225

# Write output data frames for comparison with heat maps
write.table(rep1, 
            file="spike.correlation.1-85.txt",
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

write.table(rep2, 
            file="spike.correlation.86-170.txt",
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

# Convert to matrices for heatmap function
rep1.matrix <- data.matrix(rep1)
rep2.matrix <- data.matrix(rep2)

# Create heat maps
rep1_heatmap <- heatmap(rep1.matrix, Rowv=NA, Colv=NA, col = heat.colors(256))
rep2_heatmap <- heatmap(rep2.matrix, Rowv=NA, Colv=NA, col = heat.colors(256))
