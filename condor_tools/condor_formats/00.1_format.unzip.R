#   Use readLines() to read in a list of file names, which you feed to the format.for.condor()
#       function below.
#
#   Modify the values of "formatted.vector" below appropriately
#
#   After the script is run, paste the output into an applicable condor.submit file

arguments <- commandArgs(trailingOnly=TRUE);
fastq.dir <- arguments[1];      # directory of raw files in fastq format (gzipped)
                                    # /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC/
				    # fastqs_from_core/fastqs/

# Initialize vector and get files
    formatted.vector <- NULL;
    files <- list.files(fastq.dir)

for (i in 1:length(files))   {
   
    curr.file <- files[i]
    index <- gsub(".*_S|_R.*", '', curr.file)
    read.pair <- gsub(".*_R|_001.fastq.gz", '', curr.file)
    if (read.pair == 2){ direction <- "_rev" } else { direction <- "_fwd" }

   formatted.vector[i] <- paste(
        "output=$(log_dir)/stdout_unzip.", index, direction, ".out\n",
        "error=$(log_dir)/stderr_unzip.", index, direction,".out\n",
        "log=$(log_dir)/unzip.", index, direction, ".log\n",
        "arguments=$(script_dir)/unzip.sh ", "$(data_dir)/", curr.file,
        "\nqueue 1\n",
        sep=""); 
}   #   while

output.file.name <- paste("00.1_formatted.unzip.txt", sep="");
write.table(formatted.vector,
            file=output.file.name,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE);
 


