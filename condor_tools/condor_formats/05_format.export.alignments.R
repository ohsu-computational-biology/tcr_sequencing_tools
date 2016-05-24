# Get command line arguments
arguments <- commandArgs(trailingOnly=TRUE);
list.of.files <- arguments[1];      # /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC/mixcr/alignments/

# List files from directory and sort them.
list.of.files <- list.files(list.of.files);
list.of.files <- list.of.files[order(as.numeric(gsub(".*_S|_alignment.*", '', list.of.files)))]

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
        "output=$(log_dir)stdout_mixcr_export_", index, ".out\n",
        "error=$(log_dir)stderr_mixcr_export_", index, ".out\n",
        "log=$(log_dir)mixcr_export_align_", index, ".log\n",
        "arguments=-Xmx15g -jar $(script_dir) ",
        "exportAlignments ",
	"-f ",
        "-vHit ",
        "-jHit ",
	"-dHit ",
	"-vAlignment ",
	"-jAlignment ",
	"-dAlignment ",
	"-sequence ",
	"-readId ",
	"-descrR1 ",
	"-cloneId $(data_dir)/indexes/S", index, "_index.txt ",
        "$(data_dir)align/", list.of.files[i], " ",   #   input
        "$(data_dir)export_align/", output.file.names[i], "_exported.txt", #  ouput
        "\nqueue 1\n",
        sep=""); 
}   #   while

write.table(formatted.vector,
            file="03_formatted_export_align.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE);
