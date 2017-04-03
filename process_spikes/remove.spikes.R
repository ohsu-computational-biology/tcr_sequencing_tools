##  This script removes spiked records from a fastq file

#   Inputs:
#   1.  Two text files, each of which should contain a list of records to be removed.
#       Each line in the file should include a read ID from a fastq file.  For example, 
#           a file might look like this:
#                   
#               M01416:85:000000000-AFU8J:1:1101:13087:1658 2:N:0:6
#               M01416:85:000000000-AFU8J:1:1101:16788:1694 2:N:0:6
#               M01416:85:000000000-AFU8J:1:1101:13889:1720 2:N:0:6
#               M01416:85:000000000-AFU8J:1:1101:18101:1723 2:N:0:6
#
#   2.  A fastq file
#
#   Outputs:
#   1.  A fastq file, with the specified reads removed
#   2.  A text file reporting the work done by this script

#   load depdencies
#.libPaths("/mnt/lustre1/CompBio/lib/R/library")
library(ShortRead);
library(stringr);
library(data.table);

###   TODO: handle warnings more better.  At least read.delim() is generating warnings
options(warn=-1);   

### Arguments
arguments <- commandArgs(trailingOnly=TRUE);
input.fastq <- arguments[1];
forward.reads.to.remove <- arguments[2];
out.remove.dir <- arguments[3];
out.spike.dir <- arguments[4];

###   Read in fastq file
fastq.records <- readFastq(input.fastq);
num.fastqs <- length(fastq.records);
cat(num.fastqs, " fastq reads to process\n", sep="");


###   create variables to store read ids to remove
forward.ids.to.remove <- character();

    
###   Some samples returned an empty list of IDs to remove, raising an 
###	exception in the original filtering code.  We added some checks here
###	to address this possibility.  TODO:  refine this

### Get file size
file.size.forward <- file.size(forward.reads.to.remove);

### If reads in file
if(file.size.forward > 0)	{

    ##   Read in lists of fastqs to be removed and consolidate into one variable
    forward.ids.to.remove <- fread(forward.reads.to.remove, sep = '\t')$Reads

    ##   Add alternate version in case header is different
    alt.forward.ids.to.remove <- str_replace(forward.ids.to.remove, " 1", " 2");
    forward.ids.to.remove <- c(forward.ids.to.remove, alt.forward.ids.to.remove);
}   #   fi

# Assign to new variable name
ids.to.remove <- forward.ids.to.remove

###   This conditional is another part of error-checking, for the case when
###	there are no reads to be removed (e.g. when both .reads.to.remove
###	files are empty)

if(length(ids.to.remove) > 0)	{
    ## Search ids for those we should remove
    id.filter <- srFilter(function(x)   {
	!(x@id %in% ids.to.remove)
    },   #   function(x),
    name="id.filter");

    ## Search ids for those we should keep
    id.keep <- srFilter(function(x)   {
   	(x@id %in% ids.to.remove)
    },  #   function(x)
    name="id.keep")

    ## Remove spike reads from fastq file
    output.fastq.records <- fastq.records[id.filter(fastq.records)];

    ## Extract spike reads from fastq file
    output.spike.records <- fastq.records[id.keep(fastq.records)];
}   else   {	#	none to remove, so the input fastq is the output fastq
    output.fastq.records <- fastq.records;
    output.spike.records <- NULL
}   #   else


### Update user
num.records.removed <- length(fastq.records) - length(output.fastq.records);
cat("Removing ", num.records.removed, " fastq records\n", sep="");
   
###   write the fastq out
out.fastq <- sub("[.][^.]*$", '', input.fastq)
output.file.name <- basename(paste(out.fastq, ".removed.fastq", sep=""))
output.file.name <- file.path(out.remove.dir, output.file.name)
cat("Writing output to: ", output.file.name, "\n", sep="");
writeFastq(output.fastq.records, output.file.name, compress=FALSE);

###   write the spikes out
out.spike <- sub(".assembled.fastq", '', input.fastq)
out.spike.name <- basename(paste(out.spike, ".spikes.fastq", sep=""))
out.spike.name <- file.path(out.spike.dir, out.spike.name)
cat("Writing output to: ", out.spike.name, "\n", sep = "")
writeFastq(output.spike.records, out.spike.name, compress=FALSE)

