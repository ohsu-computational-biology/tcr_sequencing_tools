#   Script for removing records with stop codons and frameshifts



remove.records <- function(input.file)    {

    #   Specify characters to search for here.
    #       Note that since "grep()" is used to identify records containing these
    #       special characters, escapes may be required (e.g. "\\*" instead of "*").
    special.characters <- c("\\*", "_");

    #   read in file
    clonotypes <- read.delim(input.file, stringsAsFactors=FALSE);
    #   identify records containing the special character
    for(i in 1:length(special.characters))  {
        aas <- clonotypes$cdr3aa;
        indices.to.remove <- integer();
        indices.to.remove <- grep(special.characters[i], aas);
        #   remove records
        if(length(indices.to.remove > 0))   {
            clonotypes <- clonotypes[-indices.to.remove,];
        }   #   fi
    }   #   for i 
    #   write output to file
    output.file.name <- paste(input.file, "_no_frameshifts_or_stop_codons.txt", sep="");

    write.table(clonotypes, 
        file=output.file.name, 
        row.names=FALSE);

}   #   remove.records()


