

#   use readLines() to read in a list of file names

format.for.condor <- function(list.of.files) {

    formatted.vector <- NULL;

    #   Assumes that there are three files per sample:
    #   
    #       S1.fastq
    #       S1.forward.reads.to.remove.txt
    #       S1.reverse.reads.to.remove.txt
  
    i <- 1; 
    j <- 1;
    loop.control <- (length(list.of.files) + 1);

    while (i < loop.control)  {

       formatted.vector[j] <- paste(
            "output=$(script_dir)logs/stdout_remove_spikes_parallel_", i, ".out\n",
            "error=$(script_dir)logs/stderr_remove_spikes_parallel_", i, ".out\n",
            "log=$(script_dir)logs/remove_spikes_parallel_", i, ".log\n",
            "arguments=$(script_dir)tcr_sequencing_tools/optimized.remove.spikes.condor.R ",
            "$(data_dir)", list.of.files[i], " ",
            "$(data_dir)", list.of.files[i + 1], " ", 
            "$(data_dir)", list.of.files[i + 2],
            "\nqueue 1\n",
            sep=""); 
       i <- i + 3;
       j <- j + 1;
    }   #   while

    write.table(formatted.vector,
                file="formatted.txt",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE);

}   #   format.for.condor()


