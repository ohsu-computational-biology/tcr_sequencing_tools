#!/usr/bin/env Rscript

#   This function calculates a "global" scaling factor
#   The factor is calculated using counts of EACH spike for ALL samples
#   The calculated scaling factor should be applied to each sample's spike count

suppressMessages(library(data.table))
suppressMessages(library(optparse))

####################
### COMMAND LINE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

optlist <- list(
  make_option(
    c("-i", "--inputFiles"),
    type = "character",
    help = "25-bp spike count files."
  ),
  make_option(
    c("-r", "--reference"),
    type = "character",
    help = "Barcode reference file with sequence and V/J identities."
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    help = "Single-column file with scaling factor for each V/J combination."
  )
)

### Parse Command Line
p <- OptionParser(usage = "%prog -i inputFiles -r reference -o output",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

### Get command-line arguments
inputs <- args$inputFiles
ref <- args$reference
output <- args$output

#########################
### INPUT AND WRANGLE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#########################

### Separate multiple input files into a list of individual files

files <- unlist(strsplit(inputs, ','))

### Read in first file to construct matrix that will be populated with count data.
### There should be one row for each spike and one column for each sample in the batch.

temp.file <- files[1]
temp <- fread(temp.file);

num.rows <- temp[,.N]
num.cols <- length(files)

counts.and.samples <- matrix(nrow = num.rows, ncol = num.cols);

### Read in each of the spike count files and append count data to counts.and.samples matrix.
for(i in 1:length(files))   {

    ## Read in file
    curr.file <- files[i]
    curr.spike.counts <- fread(curr.file)

    ## Add counts to matrix
    counts.and.samples[,i] <- curr.spike.counts$spike.count
} # for

### Add one to all cells to make sure to avoid having zeroes
counts.and.samples <- apply(counts.and.samples, c(1,2), function(x) x + 1)

############
### CALC ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############

### Calculate mean number of spikes found, across all samples.
sum.across.samples <- apply(counts.and.samples, 1, sum);
spike.count.mean <- sum(sum.across.samples) / num.rows;

### Use the mean to calculate the scaling factor.
scaling.factor <- sum.across.samples / spike.count.mean;
scaling.factor <- 1 / scaling.factor;

### Have to now add entry for V122 because we use the same spike/scaling factor for V121 and V122 #######
### Read in reference spike table
ref_spikes <- fread(ref)

### Combine V and J columns with spike counts
temp.vj <- as.data.frame(cbind(ref_spikes$V, ref_spikes$J, scaling.factor), stringsAsFactors = F)
### Extract the V121/V122 spike
v.122 <- temp.vj[temp.vj$V1 == "V12-1-2-",]
### rename V12-1-2 to V121 and V122
temp.vj$V1 <- gsub("V12-1-2-", "V12-1-", temp.vj$V1)
v.122$V1 <- gsub("V12-1-2-", "V12-2-", v.122$V1)

### Add V122 to the other scaling factors
final.vj <- as.data.frame(rbind(temp.vj, v.122), stringsAsFactors = F)

##########################
### WRANGLE AND OUTPUT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########################

### Update column names
colnames(final.vj) <- c("V", "J", "naive")

###   Remove trailing dashes
final.vj$V <- gsub("-", "", final.vj$V)

### Write output.
write.table(final.vj,
            file=output,
            quote=FALSE,
            row.names=FALSE)

