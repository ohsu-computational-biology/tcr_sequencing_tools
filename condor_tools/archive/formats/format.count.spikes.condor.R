

#   use readLines() to read in a list of file names

format.for.condor <- function(list.of.files) {

    formatted.vector <- NULL;
  
    for (i in 1:length(list.of.files))   {

       formatted.vector[i] <- paste(
            "output=$(script_dir)/logs/stdout_count_spikes_parallel_", i, ".out\n",
            "error=$(script_dir)/logs/stderr_count_spikes_parallel_", i, ".out\n",
            "log=$(script_dir)/logs/count_spikes_parallel_", i, ".log\n",
            "arguments=$(script_dir)count.spikes.wrapper.9bp.fwd.R ",
            "$(data_dir)", list.of.files[i], " ",
            "\nqueue 1\n",
            sep=""); 
    }   #   while

    write.table(formatted.vector,
                file="formatted.txt",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE);

}   #   format.for.condor()


