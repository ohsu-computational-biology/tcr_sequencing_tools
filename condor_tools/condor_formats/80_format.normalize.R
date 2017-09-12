# Get command line args
arguments <- commandArgs(trailingOnly=TRUE);
list.of.clone.files <- arguments[1];		# /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC/
		       				# normalization/decontam/
list.of.count.files <- arguments[2];		# /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC/
                                        # normalizaiton/counts/
out.dir <- arguments[3]

if (is.na(out.dir)){out.dir <- "../submits/"}

# List files from directory and sort them.
list.of.clone.files <- list.files(list.of.clone.files);
list.of.clone.files <- list.of.clone.files[order(as.numeric(gsub(".*_S|_alignment_.*", '', list.of.clone.files)))]
list.of.count.files <- list.files(list.of.count.files);
list.of.count.files <- list.of.count.files[order(as.numeric(gsub(".*_S|\\..*", '', list.of.count.files)))]

# Initialize vectors
formatted.vector <- paste("#!/bin/sh\n",
                          'getenv="True"',
                          "script_dir=$ENV(tool)/normalize/",
                          "log_dir=$ENV(data)/condor_logs/normalization/",
                          "normalization_dir=$ENV(data)/normalization/",
                          "# Program", "executable=/usr/bin/Rscript\n",
                          "# Cores", "request_cpus = 1",
                          "# Memory", "request_memory = 4 GB\n", "# Arguments\n", sep = '\n')


#   TODO - verify paralellity of input files
  
for (i in 1:length(list.of.clone.files))   {

   curr.clone <- list.of.clone.files[i]
   clone.index <- gsub(".*_S|_alignment.*", '', curr.clone)


   curr.count <- list.of.count.files[i]
   count.index <- gsub(".*_S|\\..*", '', curr.count)

   if (clone.index == count.index){

       index <- clone.index

       formatted.vector[i+1] <- paste(
           "output=$(log_dir)stdout_normalize_clonotypes_", index, ".out\n",
           "error=$(log_dir)stderr_normalize_clonotypes_", index, ".out\n",
           "log=$(log_dir)condor_normalize_clonotypes_", index, ".log\n",
           "arguments=$(script_dir)normalize.clonotype.counts.condor.R ",
           "$(normalization_dir)decontam/", list.of.clone.files[i], " ",   #   clone input
           "$(normalization_dir)counts/", list.of.count.files[i], " ",    #   count input
           "$(normalization_dir)normalized_clones/ ", #  ouput directory
           "$(normalization_dir)scaling_factor.txt ",
           "$(ref_dir)/nb.scaling.factors.txt",
           "\nqueue 1\n",
           sep="");
   } # if
}   #   for


out.file <- file.path(out.dir, "80_normalize.submit")
 write.table(formatted.vector,
             file=out.file,
             row.names = FALSE,
             col.names = FALSE,
             quote = FALSE);
