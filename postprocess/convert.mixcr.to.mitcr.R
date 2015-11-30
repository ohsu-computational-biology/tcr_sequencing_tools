#   Converts an MiXCR file to an MiTCR file
#   Useful for using Jacob's normalization tool, which was designed to
#       work with MiTCR output 

expected.mitcr.colnames <- c("Seq. Count",
                            "Percent",
                            "N Sequence",
                            "AA Sequence",
                            "V segments",
                            "J segments");
relevant.mixcr.colnames <- c("Clone.count",
                            "Clone.fraction",
                            "N..Seq..CDR3",
                            "AA..seq..CDR3",
                            "All.V.hits",
                            "All.J.hits");

convert.mixcr.to.mitcr <- function(input.file)  {

    #   Read in MiXCR file
    MiXCR.data.frame <- read.csv(input.file, 
                                sep="\t", 
                                stringsAsFactors=FALSE);
    #   Get column names in MiXCR file
    MiXCR.colnames <- colnames(MiXCR.data.frame);

    #   For all relevant columns
    for(i in 1:length(expected.mitcr.colnames)) {

        #   find the location of the appropriate column name in the MiXCR column names
        index <- which(relevant.mixcr.colnames[i] == MiXCR.colnames);
        
        #   error-check
        if(length(index) == 0)
            warning("No column in MiXCR input found for MiTCR column:  ", expected.mitcr.colnames[i]);
        #   replace the column name with the associated MiTCR column name
        MiXCR.colnames[index] <- expected.mitcr.colnames[i];

    }   #   for i

    #   Replace the old column names with new column names
    colnames(MiXCR.data.frame) <- MiXCR.colnames;

    #   Write out the modified file
    output.file <- paste(input.file, ".converted.to.mitcr.csv", sep="");
    write.table(MiXCR.data.frame, 
                output.file,
                quote=FALSE,
                sep=",",
                row.names=FALSE);

}
