# Get command line arguments.
arguments <- commandArgs(trailingOnly=TRUE);
list.of.files <- arguments[1];      # /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarinig/DNAXXXXLC/mixcr/alignments/

# Get files from directory and sort them
list.of.files <- list.files(list.of.files);
list.of.files <- list.of.files[order(as.numeric(gsub(".*_S|_alignment.vdjca", '', list.of.files)))]

# Initialize vector
formatted.vector <- paste("#!/bin/sh\n",
                          'getenv="True"',
                          "script_dir=$ENV(tool)/mixcr/mixcr-2.0/mixcr.jar",
                          "data_dir=$ENV(data)/mixcr/",
                          "input_dir=$(data_dir)/alignments",
                          "index_dir=$(data_dir)indexes",
                          "report_dir=$(data_dir)/reports/assemble",
                          "output_dir=$(data_dir)/assemblies",
                          "log_dir=$ENV(data)/condor_logs/mixcr/assemble/",
                          "# Program", "executable=/usr/bin/java/\n",
                          "# Cores", "request_cpus = 4\n",
                          "# Memory", "request_memory = 16 GB\n", "# Arguments\n", sep = '\n')


# Format vector
for (i in 1:length(list.of.files))   {
   curr.file <- list.of.files[i]
   index <- gsub(".*_S|_alignment.*", '', curr.file)
   output.file.name <- sub("[.][^.]*$", "", curr.file)
   formatted.vector[i] <- paste(
        "output=$(log_dir)stdout_mixcr_assemble_", index, ".out\n",
        "error=$(log_dir)stderr_mixcr_assemble_", index, ".out\n",
        "log=$(log_dir)mixcr_assemble_", index, ".log\n",
        "arguments=-Xmx10g -jar $(script_dir) ",
        "assemble ",
	"-f ",
	"--index $(index_dir)/S", index, "_index.txt ",
        "--report $(report_dir)/S", index, "_assemble_report.txt ",   #   report
        "--threads 4 ",
        "$(input_dir)", list.of.files[i], " ",   #   input
        "$(output_dir)", output.file.name, "_clones.clns", #  ouput
        "\nqueue 1\n",
        sep=""); 
}   #   while


out.file <- file.path(out.dir, "40_assemble.submit")
write.table(formatted.vector,
            file=out.file,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE);
