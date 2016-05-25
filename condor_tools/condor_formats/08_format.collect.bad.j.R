#   Use readLines() to read in a list of file names, which you feed to the format.for.condor()
#       function below.
#
#   Modify the values of "formatted.vector" below appropriately
#
#   After the script is run, paste the output into an applicable condor.submit file

arguments <- commandArgs(trailingOnly=TRUE);
list.of.files <- arguments[1];      # directory of despiked files in fastq format
                                    # /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC/
				    # mixcr/despiked_fasqs
files.to.extract <- arguments[2]    # directory of reads collected by mixcr_qc.R that failed due to no J alignment
		    		    # /base/mixcr/QC/align



    formatted.vector <- NULL;

    fastq.files <- list.files(list.of.files)
    fastq.sorted <- fastq.files[order(as.numeric(gsub(".*_S|\\..*", '', fastq.files)))] # sort by S## so that log numbers
    	      						    	      # correspond to sample numbers

    reads.to.extract <- list.files(files.to.extract)
    reads.to.extract <- reads.to.extract[order(as.numeric(gsub("S|_align.*", '', reads.to.extract)))]
  
    for (i in 1:length(fastq.sorted))   {
      
       # Extract index string from file, in case we delete a file due to QC...the logs will still have same number as files.
       curr.fastq <- fastq.sorted[i]
       fastq.index <- gsub(".*_S|\\..*", '', curr.fastq)

       curr.extract <- reads.to.extract[i]
       extract.index <- gsub("S|_align.*", '', curr.extract)

       if (fastq.index == extract.index) {
       
              index <- fastq.index

              formatted.vector[i] <- paste(
       	           "output=$(log_dir)/stdout_collect.bad.j.", index, ".out\n",
            	   "error=$(log_dir)/stderr_collect.bad.j.", index, ".out\n",
            	   "log=$(log_dir)/collect.bad.j.", index, ".log\n",
            	   "arguments=$(script_dir)/collect.bad.j.reads.R ", "$(fastq_dir)/", curr.fastq, " $(reads_dir)/", curr.extract, 
            	   " $(out_dir)/ ",
            	   "\nqueue 1\n",
            	   sep=""); 
	} else {
	    stop("Fastq index and reads to extract index do not match.")
	} # if
    }   #   for

    output.file.name <- "08_formatted.collect.bad.j.txt"
    write.table(formatted.vector,
                file=output.file.name,
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE);




