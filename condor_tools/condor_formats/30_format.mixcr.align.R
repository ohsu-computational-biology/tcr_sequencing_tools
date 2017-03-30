
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

  
for (i in 1:length(list.of.files))   {
   curr.file <- list.of.files[i]
   index <- gsub(".*_S|\\..*", '', curr.file)
   output.file.name <- gsub("[.a-z]", '', curr.file)
   formatted.vector[i] <- paste(
        "output=$(log_dir)stdout_mixcr_align_", index, ".out\n",
        "error=$(log_dir)stderr_mixcr_align_", index, ".out\n",
        "log=$(log_dir)mixcr_align_", index, ".log\n",
        "arguments=-Xmx15g -jar $(script_dir) ",
        "align ",
	"-f ",
        "-c TRB ",
        "--species mmu ",
	"--save-description ",
	"--save-reads ",
	"-v ",
        "--report $(report_dir)/S", index, "_align_report.txt ",   #   report
        "$(input_dir)/", list.of.files[i], " ",   #   input
        "$(out_dir)/", output.file.name, "_alignment.vdjca", #  ouput
        "\nqueue 1\n",
        sep=""); 
}   #   for

write.table(formatted.vector,
            file="03_formatted_align.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE);
