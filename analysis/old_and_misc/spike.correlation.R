########################################################################
### Script to calculate correlation of spike counts between samples ###
########################################################################

# Works if I run it in Rstudio, but not when run from command line...


# Takes a directory filled with 25 bp count files and generates a heat map based on:
#       Correlation coefficient (r) - strength and direction of linear relationship
#       Coefficient of Determination (r^2) - proportion of variance in one variable that is explained by the other

###
### Get command line arguments
###

arguments <- commandArgs(trailingOnly = TRUE);
count.directory <- arguments[1]       # Should be DNAXXXLC/spike_counts/25bp/counts/
offset <- numeric(arguments[2])                # If 170 samples, 1-85 are rep1 and 86-170 are rep2, offset = 85
comparison <- arguments[3]                  # r for correlation coefficient and r2 for coefficient of determination
type <- arguments[4]


###
### Get the data
###
# Extract list of files from directory
count.files <- list.files(count.directory);

# Order them numerically
count.files <- count.files[order(as.numeric(gsub(".*_S|\\.assembled.*", '', count.files)))]

###
### Function for replicates
###

replicates <- function(count.files, offset, comparison) {
  
  # Create data frames to hold comparison values. Should have number of rows 
  # and columns equal to number of replicates.
  
  #rep.length <- strtoi((length(count.files) - offset))
#  replicate.1 <- data.frame(matrix(ncol = rep.length, nrow = rep.length))
#      colnames(replicate.1) <- c(1:rep.length)
#      rownames(replicate.1) <- c(1:rep.length)
  replicate.1 <- data.frame(matrix(ncol = 85, nrow = 85))
    colnames(replicate.1) <- c(1:85)
    rownames(replicate.1) <- c(1:85)
  
#  replicate.2 <- data.frame(matrix(ncol = rep.length, nrow = rep.length))
#      colnames(replicate.2) <- c((rep.length+1):length(count.files))
#      rownames(replicate.2) <- c((rep.length+1):length(count.files))
    replicate.2 <- data.frame(matrix(ncol = 85, nrow = 85))
    colnames(replicate.2) <- c(86:170)
    rownames(replicate.2) <- c(86:170)
    
  # Populate data frames
  # Replicate 1
#  for (i in 1:rep.length){
  for (i in 1:85){
    curr.count <- count.files[i]
    curr.count.data <- read.delim(file.path(count.directory, curr.count), sep = ",",
                                  check.names = FALSE, stringsAsFactors = FALSE);
    curr.spikes <- curr.count.data$spike.count
#    for (j in 1:rep.length){
    for (j in 1:85){
      curr.compare.count <- count.files[j]
      curr.compare.count.data <- read.delim(file.path(count.directory, curr.compare.count), sep = ",",
                                            check.names = FALSE, stringsAsFactors = FALSE);
      curr.compare.spikes <- curr.compare.count.data$spike.count
      if (comparison == "r"){
        replicate.1[i,j] <- cor(curr.spikes, curr.compare.spikes)
      } else if (comparison == "r2") {
        replicate.1[i,j] <- summary(lm(curr.spikes ~ curr.compare.spikes))$r.squared
      } else {
        print("Error: Incorrect 'comparison' in command args. Please enter 'r' or 'r2'.")
      } # fi
    } # for j
  } # for i
  
  # Replicate 2
#  for (i in 1:rep.length){
  for (i in 1:85){
    curr.count <- count.files[i+offset]
    curr.count.data <- read.delim(file.path(count.directory, curr.count), sep = ",",
                                  check.names = FALSE, stringsAsFactors = FALSE);
    curr.spikes <- curr.count.data$spike.count
#    for (j in 1:rep.length){
    for (j in 1:85){
      curr.compare.count <- count.files[j+offset]
      curr.compare.count.data <- read.delim(file.path(count.directory, curr.compare.count), sep = ",",
                                            check.names = FALSE, stringsAsFactors = FALSE);
      curr.compare.spikes <- curr.compare.count.data$spike.count
      if (comparison == "r"){
        replicate.2[i,j] <- cor(curr.spikes, curr.compare.spikes)
      } else if (comparison == "r2") {
        replicate.2[i,j] <- summary(lm(curr.spikes ~ curr.compare.spikes))$r.squared
      } else {
        print("Error: Incorrect 'comparison' in command args. Please enter 'r' or 'r2'.")
      } # fi
    } # for j
  } # for i
      
  # Calculate average comparison value
#  row.sum.1 <- length(1:rep.length)
#  row.sum.2 <- length(1:rep.length)
  row.sum.1 <- length(1:85)
  row.sum.2 <- length(1:85)
#  for (i in 1:rep.length){
  for (i in 1:85){
    # Sum across rows for each sample
    row.sum.1[i] <- sum(replicate.1[i,])
    row.sum.2[i] <- sum(replicate.2[i,])
  } # for i
  # Aggregate by row sum and divide by total number of comparisons
#  total.comp <- rep.length * rep.length
  total.comp <- 85 * 85
  total.sum.1 <- sum(row.sum.1)
  total.sum.2 <- sum(row.sum.2)
  avg.comparison.1 <- total.sum.1 / total.comp
  avg.comparison.2 <- total.sum.2 / total.comp
  
  # Print to stdout
  if (comparison == "r"){
    print(c("Average correlation coefficient (r) for replicate 1: ", round(avg.comparison.1, 4)))
    print(c("Average correlation coefficient (r) for replicate 2: ", round(avg.comparison.2, 4)))
  } else {
    print(c("Average coefficient of determination (r^2) for replicate 1: ", round(avg.comparison.1, 4)))
    print(c("Average coefficient of determination (r^2) for replicate 2: ", round(avg.comparison.2, 4)))
  } # fi
  
  # Write output tables
  write.table(replicate.1, 
              file=paste("spike.correlation.1-", "85", ".txt", sep=''),
              quote=FALSE,
              sep="\t",
              row.names=FALSE)
  
  write.table(replicate.2, 
              file=paste("spike.correlation.", "86", "-", length(count.files), ".txt", sep = ""),
              quote=FALSE,
              sep="\t",
              row.names=FALSE)
  
  # Convert to matrices for heatmap function
  rep1.matrix <- data.matrix(replicate.1)
  rep2.matrix <- data.matrix(replicate.2)
  
  # Create heat maps
  pdf(file=file.path(count.directory, "heatmaps.pdf"))
  rep1_heatmap <- heatmap(rep1.matrix, Rowv=NA, Colv=NA, col = heat.colors(256), main =
                            "Correlations of Spike Counts Samples 1-85")
  rep2_heatmap <- heatmap(rep2.matrix, Rowv=NA, Colv=NA, col = heat.colors(256), main =
                            "Correlations of Spike Counts Samples 86-170")
  dev.off()

  
} # replicates()

no.replicates <- function(count.files, offset, comparison) {
  
  # Create data frames to hold comparison values. Should have number of rows 
  # and columns equal to number of replicates.
  
  #rep.length <- strtoi((length(count.files) - offset))
  #  replicate.1 <- data.frame(matrix(ncol = rep.length, nrow = rep.length))
  #      colnames(replicate.1) <- c(1:rep.length)
  #      rownames(replicate.1) <- c(1:rep.length)
  replicate.1 <- data.frame(matrix(ncol = 85, nrow = 85))
  colnames(replicate.1) <- c(1:85)
  rownames(replicate.1) <- c(1:85)
  
  # Populate data frames
  # Replicate 1
  #  for (i in 1:rep.length){
  for (i in 1:length(count.files)){
    curr.count <- count.files[i]
    curr.count.data <- read.delim(file.path(count.directory, curr.count), sep = ",",
                                  check.names = FALSE, stringsAsFactors = FALSE);
    curr.spikes <- curr.count.data$spike.count
    #    for (j in 1:rep.length){
    for (j in 1:length(count.files)){
      curr.compare.count <- count.files[j]
      curr.compare.count.data <- read.delim(file.path(count.directory, curr.compare.count), sep = ",",
                                            check.names = FALSE, stringsAsFactors = FALSE);
      curr.compare.spikes <- curr.compare.count.data$spike.count
      if (comparison == "r"){
        replicate.1[i,j] <- cor(curr.spikes, curr.compare.spikes)
      } else if (comparison == "r2") {
        replicate.1[i,j] <- summary(lm(curr.spikes ~ curr.compare.spikes))$r.squared
      } else {
        print("Error: Incorrect 'comparison' in command args. Please enter 'r' or 'r2'.")
      } # fi
    } # for j
  } # for i
  
  
  # Calculate average comparison value
  #  row.sum.1 <- length(1:rep.length)
  #  row.sum.2 <- length(1:rep.length)
  row.sum.1 <- numeric(length(count.files))
  #  for (i in 1:rep.length){
  for (i in 1:length(count.files)){
    # Sum across rows for each sample
    row.sum.1[i] <- sum(replicate.1[i,])
    } # for i
  # Aggregate by row sum and divide by total number of comparisons
  #  total.comp <- rep.length * rep.length
  total.comp <- length(count.files) * length(count.files)
  total.sum.1 <- sum(row.sum.1)
  avg.comparison.1 <- total.sum.1 / total.comp
    
  # Print to stdout
  if (comparison == "r"){
    print(c("Average correlation coefficient (r): ", round(avg.comparison.1, 4)))
    } else {
    print(c("Average coefficient of determination (r^2): ", round(avg.comparison.1, 4)))
    } # fi
  
  # Write output tables
  write.table(replicate.1, 
              file=paste("spike.correlation.1-", "85", ".txt", sep=''),
              quote=FALSE,
              sep="\t",
              row.names=FALSE)
      
  # Convert to matrices for heatmap function
  rep1.matrix <- data.matrix(replicate.1)
    
  # Create heat maps
  pdf(file.path(count.directory, "heatmap.pdf"))
  rep1_heatmap <- heatmap(rep1.matrix, Rowv=NA, Colv=NA, col = heat.colors(256), main =
                            "Correlation of Spike Counts")
  dev.off()
    
} # no.replicates()


# Execute
#if (type == "replicates")
do.call(type, list(count.files, offset, comparison))

