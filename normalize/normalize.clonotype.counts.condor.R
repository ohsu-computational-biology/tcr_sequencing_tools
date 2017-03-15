# Designed and developed by Jacob Bieker (jacob@bieker.tech)

# Script reads in a csv file containing spiked read counts from TCR sequencing, and calculates the 
# what number is needed to change the counts of each spiked read to the mean. 
# Using the spiked reads, it finds the corresponding VJ region in a MiTCR-formatted CSV file
# It then normalizes the count for each region in the MiTCR file using the multiples from the spikes
#
# Assumptions:
#   1.  A CSV file, named "<MiTCR File>xout.csv" per MiTCR file of the format ID,spike,count
#   2.  A MiTCR csv file per CSV file
#   3.  A CSV file detailing the barcode-to-VJ-region 
#   4.  Spiked reads are supposed to be present in the exact same frequency


#	Ensure that numeric outputs are not expressed in scientific notation
options(scipen=999);
require(data.table)

arguments <- commandArgs(trailingOnly = TRUE);
exported.clone.file <- arguments[1];  	# /DNAXXXXLC/normalization/clones/XXX_alignment_clones_exported.txt
spike.count.file <- arguments[2];	# /DNAXXXXLC/normalization/counts/XXX.assembled.spike.counts.25bp.txt
output.path <- arguments[3];		# /DNAXXXXLC/normalization/normalized_clones/
scaling.factor.file <- arguments[4];	# /DNAXXXXLC/normalization/scaling_factor.txt



### Read in the spiked_read counts
spiked_reads <- fread(spike.count.file)

### Add V122 to file and rename V12-1-2 to V121 and V122
add.V122 <- spiked_reads[spiked_reads$V == "V12-1-2-",]
spiked_reads$V <- gsub("V12-1-2-", "V12-1-", spiked_reads$V)
add.V122$V <- gsub("V12-1-2-", "V12-2-", add.V122$V)
add.V122$SPIKE_ID <- sprintf("DM_%d", 261:273)
spiked_reads <- rbind(spiked_reads, add.V122)

  
scaling.factor <- fread(scaling.factor.file)$V1

### Legacy - assign scaling.factor to spiked_reads
spiked_reads$multiples <- scaling.factor;


  # Opens the matching MiTCR file, if such file exists
#   Original commented out below
#  MiTCR_file_data <- read.csv(exported.clone.file, stringsAsFactors = FALSE)
                                        #   End original
### Read in MiXCR count data
count_data <- fread(exported.clone.file)

### THIS SHOULD BE OBSOLETE WITH ADDITION OF DECONTAM SCRIPT (DOES IT THERE INSTEAD)
### Leaving here as a reminder, in case something changes and we get weird behavior.
### Get rid of the TRB that is before every V and J segment name, so it can be matched later
###  MiTCR_file_data$"V segments" <- sub("TRB", "", MiTCR_file_data$"Best V hit");
###  MiTCR_file_data$"J segments" <- sub("TRB", "", MiTCR_file_data$"Best J hit");
###  MiTCR_file_data$"V segments" <- sub("\\*00", "", MiTCR_file_data$"V segments");
###  MiTCR_file_data$"J segments" <- sub("\\*00", "", MiTCR_file_data$"J segments");
  
### Remove the extra characters for the V segments in the spiked counts, so matches occur
spiked_reads$V <- gsub("-","", spiked_reads$V)

### Remove dashes from MiTCR data as well ALSO SHOULD BE OBSOLETE BY DECONTAM SCRIPT
###  MiTCR_file_data$`V segments` <- gsub("-", "", MiTCR_file_data$`V segments`)

### TODO:  are we right to set these to zero? 
count_data$"Normalized clone count" <- 0;
count_data$"Normalized clone fraction" <- 0;

### Go through every spike in the spike file
for(index in 1:nrow(spiked_reads)) {
    ## Get the current spike information:  spike ID, count, V-region, J-region, etc.
    current.spike.info <- spiked_reads[index,]
    ## Grab all the rows of the exported.clone.file that match both the V- and J- region
    ##     of the current spike count
    indices.to.modify <- which((count_data$`V segments` == current.spike.info$V &
                                count_data$`J segments` == current.spike.info$J))

    if(length(indices.to.modify) > 0)   {
        count_data[indices.to.modify,`Normalized clone count` := current.spike.info$multiples * count_data[indices.to.modify,`Clone count`]]
    }   #   fi
}   #   for index

###	It does not make sense to have fractions of counts, so round accordingly
count_data$"Normalized clone count" <- round(count_data$"Normalized clone count", digits=0);
###	Adjust the clone fraction
normalized.clone.count.sum <- sum(count_data$"Normalized clone count");
count_data$"Normalized clone fraction" <- count_data$"Normalized clone count" / normalized.clone.count.sum;
	
### TODO: do we still have this code? Makes VDJtools compatible
###  MiTCR_file_data <- postprocess.normalization.output(MiTCR_file_data);

### Output Names
output.file.name <- sub("[.][^.]*$", "", basename(exported.clone.file));
output.file.name <- paste(output.path, output.file.name, "_normalized.txt", sep="");

### Write table
write.table(count_data,
            output.file.name,
            quote = FALSE,
            row.names = FALSE,
            sep="\t");

### Update
cat("Writing output to: ", output.file.name, "\n");
