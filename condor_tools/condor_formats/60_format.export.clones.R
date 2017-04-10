# Get command line arguments
arguments <- commandArgs(trailingOnly=TRUE);
list.of.files <- arguments[1];      # /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC/mixcr/assemblies/
out.dir <- arguments[2]

if (is.na(out.dir)){out.dir <- "../submits/"}

# List files from directory and sort them.
list.of.files <- list.files(list.of.files);
list.of.files <- list.of.files[order(as.numeric(gsub(".*_S|_alignment.*", '', list.of.files)))]

# Initialize vectors
formatted.vector <- paste("#!/bin/sh\n",
                          'getenv="True"',
                          "script_dir=$ENV(tool)/mixcr/mixcr-2.0/mixcr.jar",
                          "data_dir=$ENV(data)/mixcr/",
                          "input_dir=$(data_dir)/assemblies/",
                          "index_dir=$(data_dir)indexes/",
                          "output_dir=$(data_dir)/export_clones/",
                          "log_dir=$ENV(data)/condor_logs/mixcr/export_clones/",
                          "# Program", "executable=/usr/bin/java\n",
                          "# Cores", "request_cpus = 1\n",
                          "# Memory", "request_memory = 12 GB\n", "# Arguments\n", sep = '\n')

# Format vector
for (i in 1:length(list.of.files))   {
   curr.file <- list.of.files[i]
   index <- gsub(".*_S|_alignment.*", '', curr.file)
   output.file.name <- sub("[.][^.]*$", "", curr.file)
   formatted.vector[i+1] <- paste(
        "output=$(log_dir)stdout_mixcr_export_clones_", index, ".out\n",
        "error=$(log_dir)stderr_mixcr_export_clones_", index, ".out\n",
        "log=$(log_dir)mixcr_export_clones_", index, ".log\n",
        "arguments=-Xmx10g -jar $(script_dir) ",
        "exportClones ",
	"-f ",
        "--filter-out-of-frames ",
        "--filter-stops ",
        "--preset full ",
	"-vAlignment ",
	"-dAlignment ",
	"-jAlignment ",
	"-vAlignments ",
	"-dAlignments ",
	"-jAlignments ",
	"-vHit ",
	"-dHit ",
	"-jHit ",
	"-readIds $(index_dir)S", index, "_index.txt ",
	"-cloneId $(index_dir)S", index, "_index.txt ",
        "$(input_dir)", list.of.files[i], " ",   #   input
        "$(output_dir)", output.file.name, "_exported.txt", #  ouput
        "\nqueue 1\n",
        sep=""); 
}   #   while


out.file <- file.path(out.dir, "60_export.clones.submit")
write.table(formatted.vector,
            file=out.file,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE);
