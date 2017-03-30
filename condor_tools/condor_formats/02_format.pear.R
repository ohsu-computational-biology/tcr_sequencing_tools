#   Use readLines() to read in a list of file names, which you feed to the format.for.condor()
#       function below.
#
#   Modify the values of "formatted.vector" below appropriately
#
#   After the script is run, paste the output into an applicable condor.submit file

arguments <- commandArgs(trailingOnly=TRUE);
list.of.files <- arguments[1];      # directory of raw files in fastq format
                                    # /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC/
				    # fastqs_from_core/fastqs/

### Initialize vector and get files
formatted.vector <- NULL;
files <- list.files(list.of.files)

### Separate into forward and reverse and order by sample number
forward <- files[gsub(".*_R", '', files) == "1_001.fastq"]
forward <- forward[order(as.numeric(gsub(".*_S|_R.*", '', forward)))]

reverse <- files[gsub(".*_R", '', files) == "2_001.fastq"]
reverse <- reverse[order(as.numeric(gsub(".*_S|_R.*", '', reverse)))]


### Begin formatting  
for (i in 1:length(forward))   {
   
   curr.forward <- forward[i]
   curr.reverse <- reverse[i]
    
   forward.index <- gsub(".*_S|_R.*", '', curr.forward)
   reverse.index <- gsub(".*_S|_R.*", '', curr.reverse)

   # Check for parallelism of files
   if (forward.index != reverse.index) {stop("Sample numbers do not match between forward and reverse.")}

   index <- forward.index

   formatted.vector[i] <- paste(
        "output=$(log_dir)/stdout_pear.", index, ".out\n",
        "error=$(log_dir)/stderr_pear.", index, ".out\n",
        "log=$(log_dir)/pear.", index, ".log\n",
        "arguments=$(script_dir)/run_pear_exacloud.pl ", "$(data_dir)/", curr.forward, " $(data_dir)/", curr.reverse,
	" -o $(out_dir)/ -f $(qc_dir)/pear_full_log.txt -s $(qc_dir)/pear_summary_log.txt",
        "\nqueue 1\n",
        sep=""); 
}   #   while

output.file.name <- paste("02_formatted.pear.txt", sep="");
write.table(formatted.vector,
            file=output.file.name,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE);
 


