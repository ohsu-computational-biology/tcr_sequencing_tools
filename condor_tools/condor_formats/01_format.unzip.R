#   Use readLines() to read in a list of file names, which you feed to the format.for.condor()
#       function below.
#
#   Modify the values of "formatted.vector" below appropriately
#
#   After the script is run, paste the output into an applicable condor.submit file

arguments <- commandArgs(trailingOnly=TRUE);
fastq.dir <- arguments[1];      # directory of raw files in fastq format (gzipped)
                                    # /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC/
				    # fastqs_from_core/fastqs/

## Initialize vector and get files
formatted.vector <- NULL;
files <- list.files(fastq.dir)

## Check parallelism
fwd_files <- files[which(gsub(".*_R|_001.fastq.gz", '', files) == "1")]
rev_files <- files[which(gsub(".*_R|_001.fastq.gz", '', files) == "2")]
if (length(fwd_files) != length(rev_files)) stop("Unequal number of fwd and rev files.")

## Get file numbers
file_nums <- as.numeric(gsub(".*_S|_R.*", '', fwd_files)); file_nums <- sort(file_nums)

## Get max sample number information
max_val <- max(file_nums)/10
max_tens <- floor(max_val)
final_ones <- as.numeric(strsplit(as.character(max_val), split = '\\.')[[1]][2])

## Special case if last file is a tens
if (is.na(final_ones)) final_ones <- 0

## Get base file information
base_name <- gsub("_S[0-9]+_.*", "_S", files[1])
suffix <- "_R[1-2]_001.fastq.gz"

for (i in 0:max_tens){
    if (i == 0){
        curr.file <- paste0(base_name, "[1-9]", suffix, collapse = '')
        cat(curr.file, "\n")
    } else if (i != max_tens){
        curr.file <- paste0(base_name, i, "[0-9]", suffix, collapse = '')
        cat(curr.file, "\n")
    } else {
        final_regex <- paste0("[0-", final_ones, "]", collapse = '')
        curr.file <- paste0(base_name, i, final_regex, suffix, collapse = '')
        cat(curr.file, "\n")
    } # fi

    formatted.vector[i+1] <- paste(
        "output=$(log_dir)/stdout_unzip.", i, ".out\n",
        "error=$(log_dir)/stderr_unzip.", i, ".out\n",
        "log=$(log_dir)/unzip.", i, ".log\n",
        "arguments=$(script_dir)/unzip.sh ", "$(data_dir)/", curr.file,
        "\nqueue 1\n",
        sep=""); 

} # for    


## for (i in 1:length(files))   {
   
##     curr.file <- files[i]
##     index <- gsub(".*_S|_R.*", '', curr.file)
##     read.pair <- gsub(".*_R|_001.fastq.gz", '', curr.file)
##     if (read.pair == 2){ direction <- "_rev" } else { direction <- "_fwd" }

##    formatted.vector[i] <- paste(
##         "output=$(log_dir)/stdout_unzip.", index, direction, ".out\n",
##         "error=$(log_dir)/stderr_unzip.", index, direction,".out\n",
##         "log=$(log_dir)/unzip.", index, direction, ".log\n",
##         "arguments=$(script_dir)/unzip.sh ", "$(data_dir)/", curr.file,
##         "\nqueue 1\n",
##         sep=""); 
## }   #   while

output.file.name <- paste("01_formatted.unzip.txt", sep="");
write.table(formatted.vector,
            file=output.file.name,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE);
 


