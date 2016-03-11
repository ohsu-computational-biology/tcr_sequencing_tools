
#   This script assumes the input files have the following format:
#
#       XXX_XX_xx.assembled.fastq.removed.fastq
#
#   (The part that really matters in the format are the "."s.)

#   
#   use readLines() to read in a list of file names

# Get command line arguments
arguments <- commandArgs(trailingOnly=TRUE);
list.of.files <- arguments[1];      # /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC/mixcr/
                                    # despiked_fastqs

# Get files and sort them
list.of.files <- list.files(list.of.files)
list.of.files <- list.of.files[order(as.numeric(gsub(".*_S|\\..*", '', list.of.files)))]

# Initialize vectors
formatted.vector <- NULL;
output.file.names <- NULL;

#   strip extensions from file names, for output file names
for(i in 1:length(list.of.files))   {
    output.file.names[i] <-  sub("[.][^.]*$", "", list.of.files[i]);
    output.file.names[i] <-  sub("[.][^.]*$", "", output.file.names[i]);
    output.file.names[i] <-  sub("[.][^.]*$", "", output.file.names[i]);
    output.file.names[i] <-  sub("[.][^.]*$", "", output.file.names[i]);
}   # for 
  
for (i in 1:length(list.of.files))   {
   formatted.vector[i] <- paste(
        "output=$(log_dir)stdout_mixcr_align_", i, ".out\n",
        "error=$(log_dir)stderr_mixcr_align_", i, ".out\n",
        "log=$(log_dir)mixcr_align_", i, ".log\n",
        "arguments=-Xmx10g -jar $(script_dir) ",
        "align ",
        "--loci TRB ",
        "--species mmu ",
        "--report $(data_dir)reports/S", i, "_align_report.txt ",   #   report
        "$(data_dir)despiked_fastqs/", list.of.files[i], " ",   #   input
        "$(data_dir)alignments/", output.file.names[i], "_alignment.vdjca", #  ouput
        "\nqueue 1\n",
        sep=""); 
}   #   for

write.table(formatted.vector,
            file="formatted_for_mixcr_align.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE);
