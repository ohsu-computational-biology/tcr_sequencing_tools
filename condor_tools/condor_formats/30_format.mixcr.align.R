
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
out.dir <- arguments[3]

# Get files and sort them
list.of.files <- list.files(list.of.files)
list.of.files <- list.of.files[order(as.numeric(gsub(".*_S|\\..*", '', list.of.files)))]

# Initialize vector
formatted.vector <- paste("#!/bin/sh\n",
                          'getenv="True"',
                          "script_dir=$ENV(tool)/mixcr/mixcr-2.0/mixcr.jar",
                          "data_dir=$ENV(data)/mixcr/",
                          "input_dir=$(data_dir)/despiked_fastqs",
                          "report_dir=$(data_dir)/reports/align",
                          "out_dir=$(data_dir)alignments",
                          "log_dir=$ENV(data)/condor_logs/mixcr/align/",
                          "# Program", "executable=/usr/bin/java/\n",
                          "# Cores", "request_cpus = 4\n",
                          "# Memory", "request_memory = 12 GB\n", "# Arguments\n", sep = '\n')


  
for (i in 1:length(list.of.files))   {
   curr.file <- list.of.files[i]
   index <- gsub(".*_S|\\..*", '', curr.file)
   output.file.name <- gsub("[.a-z]", '', curr.file)
   formatted.vector[i+1] <- paste(
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


out.file <- file.path(out.dir, "30_align.submit")
write.table(formatted.vector,
            file=out.file,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE);
