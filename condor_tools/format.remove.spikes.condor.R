

#   use readLines() to read in a list of file names

format.for.condor <- function(fastqs, reads.to.remove.forward, reads.to.remove.reverse) {

    formatted.vector <- character(length(fastqs));

    for (i in 1:length(fastqs)) {

       formatted.vector[i] <- paste(
            "output=$(script_dir)logs/stdout_remove_spikes_parallel_", i, ".out\n",
            "error=$(script_dir)logs/stderr_remove_spikes_parallel_", i, ".out\n",
            "log=$(script_dir)logs/remove_spikes_parallel_", i, ".log\n",
            "arguments=$(script_dir)tcr_sequencing_tools/optimized.remove.spikes.condor.R ",
            "$(data_dir)", fastqs[i], " ",
            "$(data_dir)", reads.to.remove.forward[i], " ", 
            "$(data_dir)", reads.to.remove.reverse[i],
            "\nqueue 1\n",
            sep=""); 
    }   #   for i

    write.table(formatted.vector,
                file="formatted.for.remove.spikes.condor.txt",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE);

}   #   format.for.condor()


