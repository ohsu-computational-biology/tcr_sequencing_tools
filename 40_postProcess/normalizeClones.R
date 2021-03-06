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

arguments <- commandArgs(trailingOnly = TRUE);
exported.clone.file <- arguments[1];  	# /DNAXXXXLC/normalization/clones/XXX_alignment_clones_exported.txt
spike.count.file <- arguments[2];	# /DNAXXXXLC/normalization/counts/XXX.assembled.spike.counts.25bp.txt
output.path <- arguments[3];		# /DNAXXXXLC/normalization/normalized_clones/
scaling.factor.file <- arguments[4];	# /DNAXXXXLC/normalization/scaling_factor.txt
nb.scaling.factor.file <- arguments[5]; # $tool/nb_norm/nb.scaling.factors.txt



### Read in the spiked_read counts
spiked_reads <- fread(spike.count.file)

### Add V122 to file and rename V12-1-2 to V121 and V122
add.V122 <- spiked_reads[spiked_reads$V == "V12-1-2-",]
spiked_reads$V <- gsub("V12-1-2-", "V12-1-", spiked_reads$V)
add.V122$V <- gsub("V12-1-2-", "V12-2-", add.V122$V)
add.V122$SPIKE_ID <- sprintf("DM_%d", 261:273)
spiked_reads <- rbind(spiked_reads, add.V122)

### Remove the extra characters for the V segments in the spiked counts, so matches occur
spiked_reads$V <- gsub("-","", spiked_reads$V)

### Read in MiXCR count data
count_data <- fread(exported.clone.file)
print(c("original row count:", count_data[,.N]))

### Remove pseudo genes
pseudo=c("V22","V8","V10","V11","V123","V18","V21","V27","V28","V25")
pseudo_reads_v <- unique(count_data$`V segments`[count_data$`V segments` %in% pseudo])
count_data <- count_data[!(count_data$`V segments` %in% pseudo),]
print(c("remove-pseudo row count:", count_data[,.N]))
print(c("removed pseudo genes test: ", pseudo_reads_v))

### Read in scaling factor file
scaling.factor_dt <- fread(scaling.factor.file)
nb.scaling.factor_dt <- fread(nb.scaling.factor.file)

print("Read in scaling factor files")

### Handle legacy code
if (ncol(scaling.factor_dt) == 1) {
    ## Notify user
    warning("The naive normalization scaling factor file only has 1 column. This file was likely ",
            "generated by an older version of calculate.scaling.factor.R. The normalization will still ",
            "proceed, but it is advised to create a new scaling factor file using the most recent version ",
            "calculate.scaling.factor.R")
    ## Extract scaling factors as vectors from their data.tables
    scaling.factor <- unlist(scaling.factor_dt[,1], use.names = F)
    nb.scaling.factor <- unlist(nb.scaling.factor_dt[,3], use.names = F)

    ## Combine spiked_reads with scaling factors
    spiked_reads$naive <- scaling.factor
    spiked_reads$nb <- nb.scaling.factor
} else {
    scaling_factors_dt <- merge(scaling.factor_dt, nb.scaling.factor_dt, by=c("V", "J"), sort = F)
    spiked_reads <- merge(spiked_reads[,c("V", "J", "spike.count")], scaling_factors_dt, by = c("V", "J"), sort = F)
    colnames(spiked_reads) <- c("V", "J", "spike.count", "naive", "nb")
} # fi

print("Done with merge")

### TODO:  are we right to set these to zero? 
count_data$"Normalized clone count" <- 0;
count_data$"Normalized clone fraction" <- 0;
count_data$"nb.clone.count" <- 0
count_data$"nb.clone.fraction" <- 0

print("Added empty columns")

### Change clone count to numeric rather than integer
countCol_v <- grep("Clone count|cloneCount", colnames(count_data), value = T)
fracCol_v <- grep("Clone fraction|cloneFraction", colnames(count_data), value = T)
changeCols_v <- c(countCol_v, fracCol_v, "Normalized clone count", "Normalized clone fraction", "nb.clone.count", "nb.clone.fraction")
print("change cols:")
print(changeCols_v)
count_data[, (changeCols_v) := lapply(.SD, as.numeric), .SDcols=changeCols_v]

### Go through every spike in the spike file
for(index in 1:nrow(spiked_reads)) {
    ## Get the current spike information:  spike ID, count, V-region, J-region, etc.
    current.spike.info <- spiked_reads[index,]
    ## Grab all the rows of the exported.clone.file that match both the V- and J- region
    ##     of the current spike count
    indices.to.modify <- which((count_data$`V segments` == current.spike.info$V &
                                count_data$`J segments` == current.spike.info$J))

    if(length(indices.to.modify) > 0)   {
        ## Create naive norm count
        count_data[indices.to.modify, `Normalized clone count` := current.spike.info$naive * count_data[indices.to.modify,get(countCol_v)]]
        ## Create nb norm count
        count_data[indices.to.modify, `nb.clone.count` := current.spike.info$nb * count_data[indices.to.modify, get(countCol_v)]]
    }   #   fi
}   #   for index

###	It does not make sense to have fractions of counts, so round accordingly
count_data$"Normalized clone count" <- round(count_data$"Normalized clone count", digits=0);
count_data$"nb.clone.count" <- round(count_data$"nb.clone.count", digits = 0)

###	Adjust the clone fraction
normalized.clone.count.sum <- sum(count_data$"Normalized clone count");
count_data$"Normalized clone fraction" <- count_data$"Normalized clone count" / normalized.clone.count.sum;

nb.clone.count.sum <- sum(count_data$"nb.clone.count")
count_data$"nb.clone.fraction" <- count_data$"nb.clone.count" / nb.clone.count.sum
	
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

warnings()
