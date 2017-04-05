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
out.dir <- arguments[3]

### Read in files and sort
bp.dir <- paste(bp, "bp", sep="");
files <- list.files(list.of.files)
sorted <- files[order(as.numeric(gsub(".*_S|\\..*", '', files)))] # sort by S## so that log numbers
                                        # correspond to sample numbers
### Initialize Vector
formatted.vector <- paste('#!/bin/sh\n',
                          'getenv="True"',
                          'script_dir=$ENV(tool)/process_spikes',
                          'data_dir=$ENV(data)/peared_fastqs/assembled',
                          paste0('out_dir=$ENV(data)/spike_counts/', bp.dir, '/'),
                          'ref_dir=$ENV(tool)/reference',
                          'log_dir=$ENV(data)/condor_logs/spike_counts/\n',
                          '# Program', 'executable=/usr/bin/Rscript\n',
                          '# Cores', 'request_cpus = 1\n',
                          '# Memory', 'request_memory = 4 GB\n', '# Arguments\n', sep = '\n')

  
for (i in 1:length(sorted))   {
    ## Extract index string from file, in case we delete a file due to QC...the logs will still have same number as files.
    curr.file <- sorted[i]
    index <- gsub(".*_S|\\..*", '', curr.file)

    formatted.vector[i+1] <- paste(
        "output=$(log_dir)/", bp.dir, "/stdout_count_spikes_parallel.", bp.dir, ".", index, ".out\n",
        "error=$(log_dir)/", bp.dir, "/stderr_count_spikes_parallel.", bp.dir, ".", index, ".out\n",
        "log=$(log_dir)/", bp.dir, "/count_spikes_parallel.", bp.dir, ".", index, ".log\n",
        "arguments=$(script_dir)/count.spikes.R ","$(data_dir)/", sorted[i], " $(ref_dir)/text_barcodesvj.txt ",
        bp, " $(out_dir)/ ",
        "\nqueue 1\n",
        sep="");
}   #   while


out.name <- paste0("10_count.spikes.", bp.dir, ".submit")
output.file.name <- file.path(out.dir, out.name);
write.table(formatted.vector,
            file=output.file.name,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE);




