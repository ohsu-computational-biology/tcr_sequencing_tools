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
		    		    # /base/mixcr/QC/j_align

d.files.to.extract <- arguments[3]  # same as above, but d_align directory



formatted.vector <- NULL;

fastq.files <- list.files(list.of.files)
fastq.sorted <- fastq.files[order(as.numeric(gsub(".*_S|\\..*", '', fastq.files)))] # sort by S## so that log numbers
    	      						    	      # correspond to sample numbers

reads.to.extract <- list.files(files.to.extract)
reads.to.extract <- reads.to.extract[order(as.numeric(gsub("S|_align.*", '', reads.to.extract)))]

d.reads.to.extract <- list.files(d.files.to.extract)
d.reads.to.extract <- d.reads.to.extract[order(as.numeric(gsub("S|_align.*", '', d.reads.to.extract)))]
  
for (i in 1:length(fastq.sorted))   {
      
   # Extract index string from file, in case we delete a file due to QC...the logs will still have same number as files.
   curr.fastq <- fastq.sorted[i]
   fastq.index <- gsub(".*_S|\\..*", '', curr.fastq)

   curr.extract <- reads.to.extract[i]
   extract.index <- gsub("S|_align.*", '', curr.extract)

   curr.d.extract <- d.reads.to.extract[i]
   d.extract.index <- gsub("S|_align.*", '', curr.d.extract)


   # Use indeces to check parallelism of files

   if (fastq.index == extract.index) { # check fastq with j ID, if doesn't match, stop and print error at end.
      if (fastq.index == d.extract.index) { # Not sure if more concise way to include both of these checks,
	     		     		      	# doing it this way to have specific error message depending on
						# which index is mismatched.
         index <- fastq.index

         formatted.vector[i] <- paste(
       	              "output=$(log_dir)/stdout_collect.bad.j.", index, ".out\n",
            	      "error=$(log_dir)/stderr_collect.bad.j.", index, ".out\n",
            	      "log=$(log_dir)/collect.bad.j.", index, ".log\n",
            	      "arguments=$(script_dir)/collect.bad.j.reads.R ", "$(fastq_dir)/", curr.fastq,
		      " $(j_dir)/", curr.extract, " $(d_dir)/", curr.d.extract, " $(j_out_dir/ ",
            	      " $(d_out_dir)/ ",
            	      "\nqueue 1\n",
            	      sep=""); 
      } else { 

         stop(paste("Fastq: ", curr.fastq, " and reads to extract: ", d.curr.extract, " do not match.", sep = ''))

      } # if fastq.index == d.extract.index

   } else {

      stop(paste("Fastq: ", curr.fastq, " and reads to extract: ", curr.extract, " do not match.", sep = ''))

   } # if fastq.index == extract.index

}   #   for

    output.file.name <- "08_formatted.collect.bad.j.txt"
    write.table(formatted.vector,
                file=output.file.name,
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE);




