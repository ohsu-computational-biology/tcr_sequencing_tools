# Get command line args
arguments <- commandArgs(trailingOnly=TRUE);
list.of.clone.files <- arguments[1];		# /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC/
		       				# normalization/decontam/
list.of.count.files <- arguments[2];		# /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC/
		       				# normalizaiton/counts/

# List files from directory and sort them.
list.of.clone.files <- list.files(list.of.clone.files);
list.of.clone.files <- list.of.clone.files[order(as.numeric(gsub(".*_S|_alignment_.*", '', list.of.clone.files)))]
list.of.count.files <- list.files(list.of.count.files);
list.of.count.files <- list.of.count.files[order(as.numeric(gsub(".*_S|\\..*", '', list.of.count.files)))]

# Initialize vectors
formatted.vector <- NULL;
output.file.names <- NULL;

#   strip extensions from file names, for output file names
for(i in 1:length(list.of.clone.files))   {
    output.file.names[i] <-  sub("[.][^.]*$", "", list.of.clone.files[i]);
}   # for 

#   TODO - verify paralellity of input files
  
for (i in 1:length(list.of.clone.files))   {

   curr.clone <- list.of.clone.files[i]
   clone.index <- gsub(".*_S|_alignment.*", '', curr.clone)


   curr.count <- list.of.count.files[i]
   count.index <- gsub(".*_S|\\..*", '', curr.count)

   if (clone.index == count.index){

       index <- clone.index

       formatted.vector[i] <- paste(
         	"output=$(log_dir)stdout_normalize_clonotypes_", index, ".out\n",
         	"error=$(log_dir)stderr_normalize_clonotypes_", index, ".out\n",
         	"log=$(log_dir)condor_normalize_clonotypes_", index, ".log\n",
         	"arguments=$(script_dir)normalize.clonotype.counts.condor.R ",
         	"$(normalization_dir)decontam/", list.of.clone.files[i], " ",   #   clone input
         	"$(normalization_dir)counts/", list.of.count.files[i], " ",    #   count input
         	"$(normalization_dir)normalized_clones/ ", #  ouput directory
         	"$(normalization_dir)scaling_factor.txt",
         	"\nqueue 1\n",
         	sep=""); 
	} # if
}   #   for

 write.table(formatted.vector,
             file="08_formatted_normalize_clonotype_counts.txt",
             row.names = FALSE,
             col.names = FALSE,
             quote = FALSE);
