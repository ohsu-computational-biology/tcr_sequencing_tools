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
Sys.getenv("vI_HOME")
if(Sys.getenv("vI_HOME")!="") {   Sys.setenv(vI_HOME="")}

library(ShortRead);
library(data.table);

###   TODO: handle warnings better.  At least read.delim() is generating warnings
options(warn=-1);   

### Arguments
arguments <- commandArgs(trailingOnly=TRUE);
input.fastq <- arguments[1];
forward.reads.to.remove <- arguments[2];
reverse.reads.to.remove <- arguments[3];
out.remove.dir <- arguments[4];
out.spike.dir <- arguments[5];
print(forward.reads.to.remove)
print(reverse.reads.to.remove)
###   Read in fastq file
fastq.records <- readFastq(input.fastq);
num.fastqs <- length(fastq.records);
cat(num.fastqs, " fastq reads to process\n", sep="");

###	create variables to store read ids to remove
forward.ids.to.remove <- character();
reverse.ids.to.remove <- character();
    
###	Some samples returned an empty list of IDs to remove, raising an
###		exception in the original filtering code.  We added some checks here
###		to address this possibility.  TODO:  refine this

### Get file size
file.size.forward <- file.size(forward.reads.to.remove);

### If reads in file
if(file.size.forward > 0)	{
    
    ##   Read in lists of fastqs to be removed and consolidate into one variable
    forward.ids.to.remove <- fread(forward.reads.to.remove, sep = '\t')$Reads
    alt.forward.ids.to.remove <- gsub(" 1", " 2", forward.ids.to.remove);
    forward.ids.to.remove <- c(forward.ids.to.remove, alt.forward.ids.to.remove);
}	#	fi
print("Forward")
print(forward.ids.to.remove[1])
print(length(forward.ids.to.remove))
### Check reverse file
file.size.reverse <- file.size(reverse.reads.to.remove);
if(file.size.reverse > 0)	{
    reverse.ids.to.remove <- fread(reverse.reads.to.remove, sep = '\t')$Reads
    alt.reverse.ids.to.remove <- gsub(" 2", " 1", reverse.ids.to.remove);
    reverse.ids.to.remove <- c(reverse.ids.to.remove, alt.reverse.ids.to.remove);
} # fi
print("Reverse")
print(reverse.ids.to.remove[1])
print(length(reverse.ids.to.remove))
### Combine forward and reverse ids
ids.to.remove <- union(forward.ids.to.remove, reverse.ids.to.remove);

###	This conditional is another part of error-checking, for the case when
###		there are no reads to be removed (e.g. when both .reads.to.remove
###		files are empty)
print("Combined")
print(ids.to.remove[1])
print(length(ids.to.remove))
if(length(ids.to.remove) > 0)	{
    ## Search ids for those we should remove
    id.filter <- srFilter(function(x)   {
        !(x@id %in% ids.to.remove)
    },   #   function(x),
    name="id.filter");

    ## Search ids for those we should keep
    id.keep <- srFilter(function(x) {
        (x@id %in% ids.to.remove)
    }, # function(x)
    name="id.keep")

    ## Remove spike reads from fastq file
    output.fastq.records <- fastq.records[id.filter(fastq.records)];

    ## Extract spike reads from fastq file
    output.spike.records <- fastq.records[id.keep(fastq.records)];
    
} else {	#	none to remove, so the input fastq is the output fastq
    output.fastq.records <- fastq.records;
    output.spike.records <- NULL
} # else

### Update user
num.records.removed <- length(fastq.records) - length(output.fastq.records);
cat("Removing ", num.records.removed, " fastq records\n", sep="");
     
###   write the fastq out
out.fastq <- gsub("_001.fastq", "", input.fastq)
output.file.name <- basename(paste(input.fastq, ".removed.fastq", sep=""));
output.file.name <- file.path(out.remove.dir, output.file.name)
cat("Writing output to: ", output.file.name, "\n", sep="");
writeFastq(output.fastq.records, output.file.name, compress=FALSE);

###   write the spikes out
out.spike <- sub("_001.fastq", '', input.fastq)
out.spike.name <- basename(paste(out.spike, ".spikes.fastq", sep=""))
out.spike.name <- file.path(out.spike.dir, out.spike.name)
cat("Writing output to: ", out.spike.name, "\n", sep = "")
writeFastq(output.spike.records, out.spike.name, compress=FALSE)


