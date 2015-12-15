

#   Use readLines() to read in a list of file names, which you feed to the format.for.condor()
#       function below.
#
#   Modify the values of "formatted.vector" below appropriately
#
#   After the script is run, paste the output into an applicable condor.submit file

format.for.condor <- function(list.of.files) {

    formatted.vector <- NULL;
  
    for (i in 1:length(list.of.files))   {

       formatted.vector[i] <- paste(
            "output=$(script_dir)/logs/stdout_count_spikes_parallel_9bp_fwd", i, ".out\n",
#            "output=$(script_dir)/logs/stdout_count_spikes_parallel_9bp_rev", i, ".out\n",
#            "output=$(script_dir)/logs/stdout_count_spikes_parallel_25bp_fwd", i, ".out\n",
            "error=$(script_dir)/logs/stderr_count_spikes_parallel_9bp_fwd", i, ".out\n",
#            "error=$(script_dir)/logs/stderr_count_spikes_parallel_9bp_rev", i, ".out\n",
#            "error=$(script_dir)/logs/stderr_count_spikes_parallel_25bp_fwd", i, ".out\n",
            "log=$(script_dir)/logs/count_spikes_parallel_9bp_fwd", i, ".log\n",
#            "log=$(script_dir)/logs/count_spikes_parallel_9bp_rev", i, ".log\n",
#            "log=$(script_dir)/logs/count_spikes_parallel_25bp_fwd", i, ".log\n",
            "arguments=$(script_dir)count.spikes.wrapper.9bp.fwd.R\n",
#            "arguments=$(script_dir)count.spikes.wrapper.9bp.rev\b",
#            "arguments=$(script_dir)count.spikes.wrapper.25bp.fwd.R\n",
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


