#   Use readLines() to read in a list of file names, which you feed to the format.for.condor()
#       function below.
#
#   Modify the values of "formatted.vector" below appropriately
#
#   After the script is run, paste the output into an applicable condor.submit file

format.for.condor <- function(list.of.files, bp="9bp", direction="fwd") {

    formatted.vector <- NULL;
    bp.and.direction <- paste(bp, ".", direction, sep="");
  
    for (i in 1:length(list.of.files))   {

       formatted.vector[i] <- paste(
            "output=$(script_dir)/logs/stdout_count_spikes_parallel.", bp.and.direction, ".", i, ".out\n",
            "error=$(script_dir)/logs/stderr_count_spikes_parallel.", bp.and.direction, ".", i, ".out\n",
            "log=$(script_dir)/logs/count_spikes_parallel.", bp.and.direction, ".", i, ".log\n",
            "arguments=$(script_dir)count.spikes.wrapper.", bp.and.direction, ".R ",
            "$(data_dir)", list.of.files[i], " ",
            "\nqueue 1\n",
            sep=""); 
    }   #   while

    output.file.name <- paste("formatted.", bp.and.direction, ".txt", sep="");
    write.table(formatted.vector,
                file=output.file.name,
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE);

}   #   format.for.condor()


