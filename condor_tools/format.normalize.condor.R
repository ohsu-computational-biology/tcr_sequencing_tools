

#   
#   use readLines() to read in lists of file names

format.for.condor <- function(list.of.clone.files, list.of.count.files) {

    formatted.vector <- NULL;
    output.file.names <- NULL;
    #   strip extensions from file names, for output file names
    for(i in 1:length(list.of.clone.files))   {
        output.file.names[i] <-  sub("[.][^.]*$", "", list.of.clone.files[i]);
    }   # for 

    #   TODO - verify paralellity of input files
  
    for (i in 1:length(list.of.clone.files))   {

       formatted.vector[i] <- paste(
            "output=$(log_dir)stdout_normalize_clonotypes_", i, ".out\n",
            "error=$(log_dir)stderr_normalize_clonotypes_", i, ".out\n",
            "log=$(log_dir)condor_normalize_clonotypes_", i, ".log\n",
            "arguments=$(script_dir)normalize.clonotype.counts.condor.R ",
            "$(normalization_dir)clones/", list.of.clone.files[i], " ",   #   clone input
            "$(normalization_dir)counts/", list.of.count.files[i], " ",    #   count input
            "$(normalization_dir)normalized_clones/ ", #  ouput directory
            "$(normalization_dir)scaling_factor.txt",
            "\nqueue 1\n",
            sep=""); 
    }   #   while

    write.table(formatted.vector,
                file="formatted_for_normalize_clonotype_counts.txt",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE);

}   #   format.for.condor()


