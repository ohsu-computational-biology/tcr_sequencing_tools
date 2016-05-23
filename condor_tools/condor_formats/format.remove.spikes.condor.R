# Get command line arguments
arguments <- commandArgs(trailingOnly=TRUE);
fastq.files <- arguments[1];        # /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC/
                                    # peared_fastqs/assembled/
to.remove.forward <- arguments[2];  # /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC/
                                    # spike_counts/9bp/reads_to_remove/

# List files and sort them.
fastqs <- list.files(fastq.files)
sorted <- fastqs[order(as.numeric(gsub(".*_S|\\..*", '', fastqs)))]

reads.to.remove.forward <- list.files(to.remove.forward)
sorted.to.remove.forward <- reads.to.remove.forward[order(as.numeric(gsub(".*_S|\\..*", '', reads.to.remove.forward)))]

# Initialize vector
formatted.vector <- character(length(fastqs));

# Format vector
for (i in 1:length(sorted)) {
   formatted.vector[i] <- paste(
        "output=$(log_dir)/stdout_remove_spikes_parallel_", i, ".out\n",
        "error=$(log_dir)/stderr_remove_spikes_parallel_", i, ".out\n",
        "log=$(log_dir)/remove_spikes_parallel_", i, ".log\n",
        "arguments=$(script_dir)tcr_sequencing_tools/process_spikes/remove.spikes.R ",
        "$(data_dir)", sorted[i], " ",
        "$(data_dir)", sorted.to.remove.forward[i], " ", 
        "$(data_dir)", "reverse.txt",
        "\nqueue 1\n",
        sep=""); 
}   #   for i

write.table(formatted.vector,
            file="formatted.for.remove.spikes.condor.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE);
