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

#scaling.factor <- calculate.scaling.factor(spike.count.dir);
#normalize.clonotype.counts(exported.clone.file, scaling.factor, sample.id);


normalize.clonotype.counts <- function(exported.clone.file, 
                                        spike.count.file, 
                                        scaling.factor, 
                                        sample.id=0)  {
  # Get the corresponding MiTCR file to go with the spiked file

  # Reads in the spiked_read counts
#   TODO:  can we make this the spike file, rather than a particular count file?
  spiked_reads <- read.csv(spike.count.file, header = TRUE, stringsAsFactors = FALSE)
  
  # Originally Jacob calculated the scaling factors in this function; the original code is
  #     reproduced below (though commented out)
  # For modularity's sake we instead calculate the scaling factor in a different function, and
  #     assign it here.  The use of "spiked_reads" is a legacy from Jacob's code 
  spiked_reads$multiples <- scaling.factor;

  # Opens the matching MiTCR file, if such file exists
#   Original commented out below
#  MiTCR_file_data <- read.csv(exported.clone.file, stringsAsFactors = FALSE)
#   End original
  MiTCR_file_data <- read.csv(exported.clone.file, stringsAsFactors = FALSE, sep="\t");
  # Get rid of the TRB that is before every V and J segment name, so it can be matched later
  # TODO:  generalize this; we assume there that there's always a trailing "*00" after the region
    #   of interest, which may not always hold, especially if MiXCR changes their output format.
    #   Similarly, we shouldn't assume there is always a leading "TRB" (unless we can validate
    #       that assumption.
  # TODO:  put some immediate error-checking in
#  MiTCR_file_data$V.segments <- gsub("^.*?V", "V", MiTCR_file_data$Best.V.hit)
#  MiTCR_file_data$J.segments <- gsub("^.*?J", "J", MiTCR_file_data$Best.J.hit)
  MiTCR_file_data$V.segments <- sub("TRB", "", MiTCR_file_data$Best.V.hit);
  MiTCR_file_data$J.segments <- sub("TRB", "", MiTCR_file_data$Best.J.hit);
  MiTCR_file_data$V.segments <- sub("\\*00", "", MiTCR_file_data$V.segments);
  MiTCR_file_data$J.segments <- sub("\\*00", "", MiTCR_file_data$J.segments);
  
  # Remove the extra characters for the V segments in the spiked counts, so matches occur
  spiked_reads$V <- gsub("-","", spiked_reads$V)

  # TODO:  are we right to set these to zero? 
  MiTCR_file_data$normalized.count <- 0; 
  MiTCR_file_data$normalized.percent <- 0; 

  # Go through every spike in the spike file
    for(index in 1:nrow(spiked_reads)) {
          # Get the current spike information:  spike ID, count, V-region, J-region, etc.
        current.spike.info <- spiked_reads[index,]
          # Grab all the rows of the exported.clone.file that match both the V- and J- region
          #     of the current spike count
        indices.to.modify <- which((MiTCR_file_data$V.segments == current.spike.info$V) & (MiTCR_file_data$J.segments == current.spike.info$J));
        if(length(indices.to.modify) > 0)   {
            MiTCR_file_data[indices.to.modify,]$normalized.count <- current.spike.info$multiples * MiTCR_file_data[indices.to.modify,]$Clone.count;
          MiTCR_file_data[indices.to.modify,]$normalized.percent <- current.spike.info$multiples * MiTCR_file_data[indices.to.modify,]$Clone.fraction;
        }   #   fi
    }   #   for index
 
    output.file.name <- paste("S", sample.id, "_exported_clones_normalized_unconverted.txt", sep="");
    write.table(MiTCR_file_data, 
                output.file.name, 
                quote = FALSE, 
                row.names = FALSE, 
                sep="\t");
    
    #   make column-names VDJTools-compatible 
    MiTCR_file_data <- postprocess.normalization.output(MiTCR_file_data);
    output.file.name <- paste("S", sample.id, "_exported_clones_normalized.txt", sep="");
    write.table(MiTCR_file_data, 
                output.file.name, 
                quote = FALSE, 
                row.names = FALSE, 
                sep="\t");
 
 
}   #   normalize.clonotype.counts()


#   TODO:  remove the dependency on column order, since this isn't guaranteed
#   TODO:  add some immediate checks to verify that the order is as we expect it to be
postprocess.normalization.output <- function(input.table)   {

    #  delete various columns that VDJTools does not expect (and possibly cannot handle)
    input.table$Clone.count <- NULL;
    input.table$Clone.fraction <- NULL;
    input.table$V.segments <- NULL;
    input.table$J.segments <- NULL;

    new.colnames <- c("Clonal sequence(s)",
                      "Clonal sequence quality(s)",
                      "All V hits",
                      "All D hits",
                      "All J hits",
                      "All C hits",  
                      "All V alignment",
                      "All D alignment",
                      "All J alignment",
                      "All C alignment",
                      "N. Seq. FR1",
                      "Min. qual. FR1",
                      "N. Seq. CDR1",
                      "Min. qual. CDR1",
                      "N. Seq. FR2",
                      "Min. qual. FR2",
                      "N. Seq. CDR2",
                      "Min. qual. CDR2",
                      "N. Seq. FR3",
                      "Min. qual. FR3",
                      "N. Seq. CDR3",
                      "Min. qual. CDR3",
                      "N. Seq. FR4",
                      "Min. qual. FR4",
                      "AA. seq. FR1",
                      "AA. seq. CDR1",
                      "AA. seq. FR2",
                      "AA. seq. CDR2",
                      "AA. seq. FR3",
                      "AA. seq. CDR3",
                      "AA. seq. FR4",
                      "Best V hit",
                      " Best J hit", 
                      "Clone count", 
                      "Clone fraction");

    #   assign new column names to data frame
    colnames(input.table) <- new.colnames;

    #   convert floating point counts to integers
    input.table$"Clone count" <- round(input.table$"Clone count");
    #   convert "Inf" values to zeroes
    indices.to.change <- which(input.table$"Clone fraction" == Inf);
    if(length(indices.to.change) > 0)   {
        input.table[indices.to.change,]$"Clone fraction" <- 0;
    }   #   fi
    indices.to.change <- which(input.table$"Clone count" == Inf);
    if(length(indices.to.change) > 0)   {
        input.table[indices.to.change,]$"Clone count" <- 0;
    }   #   fi

    return(input.table);

}   #   postprocess.normalization.output()   

#   This function calculates a "global" scaling factor
#   The factor is calculated using counts of EACH spike for ALL samples
#   The calculated scaling factor should be applied to each sample's spike count

calculate.scaling.factor <- function(spike.count.dir)  {
    
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

    return(scaling.factor);

}   #   calculate.scaling.factor()

  
  # BEGIN Jacob's original code
  #Get the mean from the last column, which is the read count
#  spiked_mean <- mean(spiked_reads[[5]])
  # Test vector holding all the multiples needed to hit the mean
#  multiples_needed <- spiked_mean/spiked_reads$spike.count
  #Puts the spiked_reads in the spiked_reads.frame for later use
#  spiked_reads$multiples <- multiples_needed
  # END Jacob's original code
