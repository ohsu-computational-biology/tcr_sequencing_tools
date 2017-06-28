
###   This function calculates a "global" scaling factor
###   The factor is calculated using counts of EACH spike for ALL samples
###   The calculated scaling factor should be applied to each sample's spike count


### Get command line arguments
library(data.table)
arguments <- commandArgs(trailingOnly=TRUE);
spike.count.dir <- arguments[1];	# DNAXXXXLC/normalization/counts/
ref.file <- arguments[2];		# /tcr_sequencing_tools/reference/text_barcodesvj.txt
    
###   Get list of files in directory
spike.count.files <- list.files(spike.count.dir);

### Confirm that files are 25bp-count files
suffixes <- sapply(spike.count.files, function(x) gsub("^[A-Z0-9]+_S[0-9]+", "", x), USE.NAMES=F)
suffixes <- unique(suffixes)
if (length(suffixes) != 1 || suffixes != ".assembled.spike.counts.25bp.txt") stop("Wrong spike input files")

###   Create matrix to hold results
###   One row for each spike; we take a peek at a spike.count.txt file to
###       find out how many spikes there are.  This aids portability, since
###       the user doesn't ned to know how many spikes there are prior to
###       running the script

temp.file <- paste(spike.count.dir, spike.count.files[1], sep="");
temp <- fread(temp.file);
num.rows <- temp[,.N]

###   One column for each sample
num.cols <- length(spike.count.files);
counts.and.samples <- matrix(nrow = num.rows, ncol = num.cols);

###   Read in each of the spike count files
for(i in 1:length(spike.count.files))   {
    ##   Read in spike.count file
    curr.file <- paste(spike.count.dir, spike.count.files[i], sep="");
    curr.spike.counts <- fread(curr.file);
    ##   Add counts to matrix
    counts.and.samples[,i] <- curr.spike.counts$spike.count;
}   #   for

## Add one to all cells to make sure to avoid having zeroes
counts.and.samples <- apply(counts.and.samples, c(1,2), function(x) x + 1)

###   Calculate mean number of spikes found, across all samples
sum.across.samples <- apply(counts.and.samples, 1, sum);
spike.count.mean <- sum(sum.across.samples) / num.rows;

###   Use the mean to calculate the scaling factor
scaling.factor <- sum.across.samples / spike.count.mean;
scaling.factor <- 1 / scaling.factor;

####### Have to now add entry for V122 because we use the same spike/scaling factor for V121 and V122 #######

###   Read in reference spike table
ref_spikes <- fread(ref.file)

###	Combine V and J columns with spike counts
temp.vj <- as.data.frame(cbind(ref_spikes$V, ref_spikes$J, scaling.factor), stringsAsFactors = F)
###   Extract the V121/V122 spike
v.122 <- temp.vj[temp.vj$V1 == "V12-1-2-",]
###   rename V12-1-2 to V121 and V122
temp.vj$V1 <- gsub("V12-1-2-", "V12-1-", temp.vj$V1)
v.122$V1 <- gsub("V12-1-2-", "V12-2-", v.122$V1)

###   Add V122 to the other scaling factors
final.vj <- as.data.frame(rbind(temp.vj, v.122), stringsAsFactors = F)

### Column names
colnames(final.vj) <- c("V", "J", "naive")

### Remove trailing dash
final.vj$V <- gsub("-", "", final.vj$V)

###   Extract just scaling factor values
#scaling.factor <- final.vj$scaling.factor

write.table(final.vj,
            file="scaling_factor.txt",
            quote=FALSE,
            row.names=FALSE);

