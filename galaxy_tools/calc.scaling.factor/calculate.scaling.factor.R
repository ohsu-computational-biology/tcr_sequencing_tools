#   This function calculates a "global" scaling factor
#   The factor is calculated using counts of EACH spike for ALL samples
#   The calculated scaling factor should be applied to each sample's spike count


#   Get command line arguments

inputs <- commandArgs(trailingOnly=TRUE)[1]
output <- commandArgs(trailingOnly=TRUE)[2]

#   Separate multiple input files into a list of individual files

files <- strsplit(inputs, ',')

#   Read in first file to construct matrix that will be populated with count data.
#   There should be one row for each spike and one column for each sample in the batch.

temp.file <- files[[1]][1]
temp <- read.csv(temp.file);

num.rows <- length(temp[[1]]);
num.cols <- length(files[[1]])

counts.and.samples <- matrix(nrow = num.rows, ncol = num.cols);


#   Read in each of the spike count files and append count data to counts.and.samples matrix.

for(i in 1:length(files[[1]]))   {

    #   Read in file
    curr.file <- files[[1]][i]
    curr.spike.counts <- read.csv(curr.file)

    #   Add counts to matrix
    counts.and.samples[,i] <- curr.spike.counts$spike.count
}   #   for

# Add one to all cells to make sure to avoid having zeroes
counts.and.samples <- apply(counts.and.samples, c(1,2), function(x) x + 1)

#   Calculate mean number of spikes found, across all samples.
#   The "apply" function calculates the sum across rows (across samples). [1 specifies row]

sum.across.samples <- apply(counts.and.samples, 1, sum);
spike.count.mean <- sum(sum.across.samples) / num.rows;

#   Use the mean to calculate the scaling factor.

scaling.factor <- sum.across.samples / spike.count.mean;
scaling.factor <- 1 / scaling.factor;

#   Write output.
write.table(scaling.factor,
            file=output,
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE);
