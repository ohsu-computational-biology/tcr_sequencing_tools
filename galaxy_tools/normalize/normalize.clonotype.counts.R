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
suppressMessages(require(data.table))


# Galaxy Argument Test
exported.clone.file <- commandArgs(trailingOnly = TRUE)[1]
spike.count.file <- commandArgs(trailingOnly = TRUE)[2]
scaling.factor.file <- commandArgs(trailingOnly = TRUE)[3]
output.file <- commandArgs(trailingOnly = TRUE)[4]


### Get spikes
spiked_reads <- fread(spike.count.file)

### Have to add V122 and rename V12-1-2 to V121
add.V122 <- spiked_reads[spiked_reads$V == "V12-1-2-",]
spiked_reads$V <- gsub("V12-1-2-", "V12-1-", spiked_reads$V)
add.V122$V <- gsub("V12-1-2-", "V12-2-", add.V122$V)
add.V122$SPIKE_ID <- sprintf("DM_%d", 261:273)
spiked_reads <- rbind(spiked_reads, add.V122)

### Read in scaling factor
scaling.factor <- as.numeric(read.table(scaling.factor.file)[,1]);

### assign it here.  The use of "spiked_reads" is a legacy from Jacob's code 
spiked_reads$multiples <- scaling.factor;

### Read in count data
count_data <- fread(exported.clone.file)
  
### Remove the extra characters for the V segments in the spiked counts, so matches occur
spiked_reads$V <- gsub("-","", spiked_reads$V)

### Initialize empty variables
count_data$"Normalized clone count" <- 0;
count_data$"Normalized clone fraction" <- 0;

### Go through every spike in the spike file
for(index in 1:nrow(spiked_reads)) {
    
    ## Get the current spike information:  spike ID, count, V-region, J-region, etc.
    current.spike.info <- spiked_reads[index,]
    ## Grab all the rows of the exported.clone.file that match both the V- and J- region
    ##     of the current spike count
    indices.to.modify <- which((count_data$`V segments` == current.spike.info$V) &
                               (count_data$`J segments` == current.spike.info$J));

    if(length(indices.to.modify) > 0)   {
        count_data[indices.to.modify, `Normalized clone count` := current.spike.info$multiples * count_data[indices.to.modify, `Clone count`]]
    }   #   fi
}   #   for index

###	It does not make sense to have fractions of counts, so round accordingly
count_data$"Normalized clone count" <- round(count_data$"Normalized clone count", digits=0);

###	Adjust the clone fraction
normalized.clone.count.sum <- sum(count_data$"Normalized clone count");
count_data$"Normalized clone fraction" <- count_data$"Normalized clone count" / normalized.clone.count.sum;

### Write output
write.table(count_data,
            file = output.file,
            quote = FALSE,
            row.names = FALSE,
            sep="\t");

 
 
