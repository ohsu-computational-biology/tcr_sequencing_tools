
#   This function calculates a "global" scaling factor
#   The factor is calculated using counts of EACH spike for ALL samples
#   The calculated scaling factor should be applied to each sample's spike count

# calculate.scaling.factor <- function(spike.count.dir)  {

# Get command line arguments
arguments <- commandArgs(trailingOnly=TRUE);
spike.count.dir <- arguments[1];	# DNAXXXXLC/normalization/counts/
ref.file <- arguments[2];		# /tcr_sequencing_tools/reference/text_barcodesvj.txt
    
    #   Get list of files in directory
    spike.count.files <- list.files(spike.count.dir);

    #   TODO:  error-check that they all have the appropriate suffix to be
    #       spike.count.txt files
    #   Create matrix to hold results
    #   One row for each spike; we take a peek at a spike.count.txt file to
    #       find out how many spikes there are.  This aids portability, since
    #       the user doesn't ned to know how many spikes there are prior to
    #       running the script
    temp.file <- paste(spike.count.dir, spike.count.files[1], sep="");
    temp <- read.csv(temp.file);
    num.rows <- length(temp[[1]]);
    #   One column for each sample
    num.cols <- length(spike.count.files);
    counts.and.samples <- matrix(nrow = num.rows, ncol = num.cols); 

    #   Read in each of the spike count files
    for(i in 1:length(spike.count.files))   {
        #   Read in spike.count file
        curr.file <- paste(spike.count.dir, spike.count.files[i], sep="");
        curr.spike.counts <- read.csv(curr.file);
        #   Add counts to matrix
        counts.and.samples[, i] <- curr.spike.counts$spike.count;
    }   #   for 

    #   Calculate mean number of spikes found, across all samples
    #   The "apply" function calculates the sum across rows (across samples)
    sum.across.samples <- apply(counts.and.samples, 1, sum);
    spike.count.mean <- sum(sum.across.samples) / num.rows;

    #   Use the mean to calculate the scaling factor
    scaling.factor <- sum.across.samples / spike.count.mean;
    scaling.factor <- 1 / scaling.factor;

####### Have to now add entry for V122 because we use the same spike/scaling factor for V121 and V122 #######

    #   Read in reference spike table
    ref_spikes <- read.table(ref.file, header = T, sep = ' ', stringsAsFactors = F)

    #	Combine V and J columns with spike counts
    temp.vj <- as.data.frame(cbind(ref_spikes$V, ref_spikes$J, scaling.factor), stringsAsFactors = F)
    #   Extract the V121/V122 spike
    v.122 <- temp.vj[temp.vj$V1 == "V12-1-2-",]
    #   rename V12-1-2 to V121 and V122
    temp.vj$V1 <- gsub("V12-1-2-", "V12-1-", temp.vj$V1)
    v.122$V1 <- gsub("V12-1-2-", "V12-2-", v.122$V1)
    #   Add V122 to the other scaling factors
    final.vj <- as.data.frame(rbind(temp.vj, v.122), stringsAsFactors = F)
    #   Extract just scaling factor values
    scaling.factor <- final.vj$scaling.factor

    write.table(scaling.factor,
                file="scaling_factor.txt",
                quote=FALSE,
                row.names=FALSE,
                col.names=FALSE);

#}   #   calculate.scaling.factor()
