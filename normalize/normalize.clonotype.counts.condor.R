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

arguments <- commandArgs(trailingOnly = TRUE);
exported.clone.file <- arguments[1];  	# /DNAXXXXLC/normalization/clones/XXX_alignment_clones_exported.txt
spike.count.file <- arguments[2];	# /DNAXXXXLC/normalization/counts/XXX.assembled.spike.counts.25bp.txt
output.path <- arguments[3];		# /DNAXXXXLC/normalization/normalized_clones/
scaling.factor.file <- arguments[4];	# /DNAXXXXLC/normalization/scaling_factor.txt



  # Reads in the spiked_read counts
#   TODO:  can we make this the spike file, rather than a particular count file?
  spiked_reads <- read.csv(spike.count.file, 
  							header = TRUE, 
							stringsAsFactors = FALSE);
  
  # Originally Jacob calculated the scaling factors in this function; the original code is
  #     reproduced below (though commented out)
  # For modularity's sake we instead calculate the scaling factor in a different function, and
  # read in scaling factor from file
  scaling.factor <- as.numeric(read.table(scaling.factor.file)[,1]);
  #     assign it here.  The use of "spiked_reads" is a legacy from Jacob's code 
  spiked_reads$multiples <- scaling.factor;

  # Opens the matching MiTCR file, if such file exists
#   Original commented out below
#  MiTCR_file_data <- read.csv(exported.clone.file, stringsAsFactors = FALSE)
#   End original
  MiTCR_file_data <- read.csv(exported.clone.file, 
  								stringsAsFactors = FALSE, 
								sep="\t",
								check.names=FALSE);
  # Get rid of the TRB that is before every V and J segment name, so it can be matched later
  # TODO:  generalize this; we assume there that there's always a trailing "*00" after the region
    #   of interest, which may not always hold, especially if MiXCR changes their output format.
    #   Similarly, we shouldn't assume there is always a leading "TRB" (unless we can validate
    #       that assumption.
  # TODO:  put some immediate error-checking in
#  MiTCR_file_data$V.segments <- gsub("^.*?V", "V", MiTCR_file_data$Best.V.hit)
#  MiTCR_file_data$J.segments <- gsub("^.*?J", "J", MiTCR_file_data$Best.J.hit)
  MiTCR_file_data$"V segments" <- sub("TRB", "", MiTCR_file_data$"Best V hit");
  MiTCR_file_data$"J segments" <- sub("TRB", "", MiTCR_file_data$"Best J hit");
  #	TODO:  do these next two lines do any work??
  MiTCR_file_data$"V segments" <- sub("\\*00", "", MiTCR_file_data$"V segments");
  MiTCR_file_data$"J segments" <- sub("\\*00", "", MiTCR_file_data$"J segments");
  
  # Remove the extra characters for the V segments in the spiked counts, so matches occur
  spiked_reads$V <- gsub("-","", spiked_reads$V)
  # Change V1212 to V121. TODO: figure out a way to make this adaptable to any situaiton.
  spiked_reads$V <- gsbu("V1212", "V121", spiked_reads$V)

  # Remove the dashes from the MiTCR data so that all V's correspond
  MiTCR_file_data$`V segments` <- gsub("-", "", MiTCR_file_data$`V segments`)

  # TODO:  are we right to set these to zero? 
  MiTCR_file_data$"Normalized clone count" <- 0; 
  MiTCR_file_data$"Normalized clone fraction" <- 0; 

  # Go through every spike in the spike file
    for(index in 1:nrow(spiked_reads)) {
          # Get the current spike information:  spike ID, count, V-region, J-region, etc.
        current.spike.info <- spiked_reads[index,]
          # Grab all the rows of the exported.clone.file that match both the V- and J- region
          #     of the current spike count
        indices.to.modify <- which((MiTCR_file_data$"V segments" == current.spike.info$V) & (MiTCR_file_data$"J segments" == current.spike.info$J));
        if(length(indices.to.modify) > 0)   {
            MiTCR_file_data[indices.to.modify,]$"Normalized clone count" <- current.spike.info$multiples * MiTCR_file_data[indices.to.modify,]$"Clone count";
		  #	This is the original code supplied by Jacob
		  #		The normalized percent (fraction) does not sum to 1 after this 
		  #		calculation, so we are commenting out and trying our own technique
#          MiTCR_file_data[indices.to.modify,]$"normalized percent" <- current.spike.info$multiples * MiTCR_file_data[indices.to.modify,]$"Clone fraction";
        }   #   fi
    }   #   for index

	#	It does not make sense to have fractions of counts, so round accordingly
	MiTCR_file_data$"Normalized clone count" <- round(MiTCR_file_data$"Normalized clone count", digits=0);
	#	Adjust the clone fraction
	normalized.clone.count.sum <- sum(MiTCR_file_data$"Normalized clone count");
	MiTCR_file_data$"Normalized clone fraction" <- MiTCR_file_data$"Normalized clone count" / normalized.clone.count.sum;
	

#   TODO:  Keep this?  Originally this output the "raw" data, not cleaning it up for use with other
#       tools.  I am not sure if there's any reason to keep this code around though, since I cannot
#       recall ever having used the "raw" data.
#   TODO:  fix file name 
#   TODO:  write regression test for this
#    output.file.name <- paste("S_exported_clones_normalized_unconverted.txt", sep="");
#    write.table(MiTCR_file_data, 
#                output.file.name, 
#                quote = FALSE, 
#                row.names = FALSE, 
#                sep="\t");
   
    #   make column-names VDJTools-compatible 
#    MiTCR_file_data <- postprocess.normalization.output(MiTCR_file_data);
	output.file.name <- sub("[.][^.]*$", "", basename(exported.clone.file));
    output.file.name <- paste(output.path, output.file.name, "_normalized.txt", sep="");
	cat("Writing output to: ", output.file.name, "\n");
    write.table(MiTCR_file_data, 
                output.file.name, 
                quote = FALSE, 
                row.names = FALSE, 
                sep="\t");
 
 


  
  # BEGIN Jacob's original code
  #Get the mean from the last column, which is the read count
#  spiked_mean <- mean(spiked_reads[[5]])
  # Test vector holding all the multiples needed to hit the mean
#  multiples_needed <- spiked_mean/spiked_reads$spike.count
  #Puts the spiked_reads in the spiked_reads.frame for later use
#  spiked_reads$multiples <- multiples_needed
  # END Jacob's original code

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

