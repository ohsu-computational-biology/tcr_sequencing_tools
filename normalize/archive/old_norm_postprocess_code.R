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
