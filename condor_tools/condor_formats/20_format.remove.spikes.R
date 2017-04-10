# Get command line arguments
arguments <- commandArgs(trailingOnly=TRUE);
fastq.files <- arguments[1];        # /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC/
                                    # peared_fastqs/assembled/
to.remove.forward <- arguments[2];  # /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC/
                                        # spike_counts/9bp/reads_to_remove/
out.dir <- arguments[3]

if (is.na(out.dir)){out.dir <- "../submits/"}

# List files and sort them.
fastqs <- list.files(fastq.files)
sorted <- fastqs[order(as.numeric(gsub(".*_S|\\..*", '', fastqs)))]

reads.to.remove.forward <- list.files(to.remove.forward)
sorted.to.remove.forward <- reads.to.remove.forward[order(as.numeric(gsub(".*_S|\\..*", '', reads.to.remove.forward)))]

# Initialize vector
formatted.vector <- paste('#!/bin/sh\n',
                          'getenv="True"',
                          'script_dir=$ENV(tool)/',
                          'data_dir=$ENV(data)/peared_fastqs/assembled/',
                          'log_dir=$ENV(data)/condor_logs/mixcr/despiked',
                          'remove_dir=$ENV(data)/spike_counts/9bp/reads_to_remove/',
                          'remove_out_dir=$ENV(data)/mixcr/despiked_fastqs/',
                          'spike_out_dir=$ENV(data)/spike_counts/9bp/spikes/\n',
                          '# Program', 'executable=/usr/bin/Rscript\n',
                          '# Cores', 'request_cpus = 1\n',
                          '# Memory', 'request_memory = 12 GB\n', '# Arguments\n', sep = '\n')
                          

# Format vector
for (i in 1:length(sorted)) {
    curr.file <- sorted[i]
    curr.remove <- sorted.to.remove.forward[i]
    
    fastq.index <- gsub(".*_S|\\..*", '', curr.file)
    reads.index <- gsub(".*_S|\\..*", '', curr.remove)
    if (fastq.index != reads.index) {stop(c("Files don't match ", i))}
    index <- fastq.index
    
   formatted.vector[i+1] <- paste(
        "output=$(log_dir)/stdout_remove_spikes_parallel_", index, ".out\n",
        "error=$(log_dir)/stderr_remove_spikes_parallel_", index, ".out\n",
        "log=$(log_dir)/remove_spikes_parallel_", index, ".log\n",
        "arguments=$(script_dir)process_spikes/remove.spikes.R ",
        "$(data_dir)", curr.file, " ",
        "$(remove_dir)", curr.remove, " $(remove_out_dir) $(spike_out_dir)", 
        "\nqueue 1\n",
        sep=""); 
}   #   for i


out.file <- file.path(out.dir, "20_remove.spikes.submit")

write.table(formatted.vector,
            file=out.file,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE);
