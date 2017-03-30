#   Use readLines() to read in a list of file names, which you feed to the format.for.condor()
#       function below.
#
#   Modify the values of "formatted.vector" below appropriately
#
#   After the script is run, paste the output into an applicable condor.submit file

arguments <- commandArgs(trailingOnly=TRUE);
list.of.files <- arguments[1];      # directory of raw files in fastq format
                                    # /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC/
				    # peared_fastqs/assembled/
bp <- arguments[2];                 # length of spike - 25 or 9
#direction <- arguments[3]           # should be FWD or REV (caps)


    formatted.vector <- NULL;
    bp.dir <- paste(bp, "bp", sep="");
#    bp.and.direction <- paste(bp, ".", direction, sep="");
    files <- list.files(list.of.files)
    sorted <- files[order(as.numeric(gsub(".*_S|\\..*", '', files)))] # sort by S## so that log numbers
    	      						    	      # correspond to sample numbers
  
    for (i in 1:length(sorted))   {
      
       # Extract index string from file, in case we delete a file due to QC...the logs will still have same number as files.
       curr.file <- sorted[i]
       index <- gsub(".*_S|\\..*", '', curr.file)

       formatted.vector[i] <- paste(
            "output=$(log_dir)/", bp.dir, "/stdout_count_spikes_parallel.", bp.dir, ".", index, ".out\n",
            "error=$(log_dir)/", bp.dir, "/stderr_count_spikes_parallel.", bp.dir, ".", index, ".out\n",
            "log=$(log_dir)/", bp.dir, "/count_spikes_parallel.", bp.dir, ".", index, ".log\n",
            "arguments=$(script_dir)/count.spikes.R ","$(data_dir)/", sorted[i], " $(ref_dir)/text_barcodesvj.txt ", 
            bp, " $(out_dir)/ ",
            "\nqueue 1\n",
            sep=""); 
    }   #   while

    output.file.name <- paste("10_formatted.", bp.dir, ".txt", sep="");
    write.table(formatted.vector,
                file=output.file.name,
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE);




