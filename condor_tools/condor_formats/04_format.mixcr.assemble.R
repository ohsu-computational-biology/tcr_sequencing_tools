# Get command line arguments.
arguments <- commandArgs(trailingOnly=TRUE);
list.of.files <- arguments[1];      # /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarinig/DNAXXXXLC/mixcr/alignments/

# Get files from directory and sort them
list.of.files <- list.files(list.of.files);
list.of.files <- list.of.files[order(as.numeric(gsub(".*_S|_alignment.vdjca", '', list.of.files)))]

# Initialize vectors
formatted.vector <- NULL;
output.file.names <- NULL;

#   strip extensions from file names, for output file names
for(i in 1:length(list.of.files))   {
    output.file.names[i] <-  sub("[.][^.]*$", "", list.of.files[i]);
}   # for 

# Format vector
for (i in 1:length(list.of.files))   {
   curr.file <- list.of.files[i]
   index <- gsub(".*_S|_alignment.*", '', curr.file)
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
        "$(output_dir)", output.file.names[i], "_clones.clns", #  ouput
        "\nqueue 1\n",
        sep=""); 
}   #   while

write.table(formatted.vector,
            file="04_formatted_assemble.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE);
