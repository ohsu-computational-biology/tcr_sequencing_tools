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
suppressMessages(source("https://bioconductor.org/biocLite.R", echo = FALSE, verbose = FALSE))
suppressMessages(library(ShortRead))
suppressMessages(library(stringr))
suppressMessages(library(data.table))

#   TODO: handle warnings more better.  At least read.delim() is generating warnings
# Setting warn to negative value ignores all warnings
options(warn=-1);   

arguments <- commandArgs(trailingOnly=TRUE);
input.fastq <- arguments[1];
reads.to.remove <- arguments[2];
output.fastq <- arguments[3]
output.spikes <- arguments[4]


#   Read in fastq file
fastq.records <- readFastq(input.fastq);
num.fastqs <- length(fastq.records);

#	create variables to store read ids to remove
ids.to.remove <- character();

#	Some samples returned an empty list of IDs to remove, raising an 
#		exception in the original filtering code.  We added some checks here
#		to address this possibility.  TODO:  refine this
file.size.to.remove <- file.size(reads.to.remove);
if(file.size.to.remove > 0)	{
#   Read in lists of fastqs to be removed and consolidate into one variable
#   TODO:  find a cleaner way to import
ids.to.remove <- unlist(fread(reads.to.remove, sep = '\t'), use.names = F)

alt.ids.to.remove <- str_replace(ids.to.remove, " 1", " 2");
ids.to.remove <- c(ids.to.remove, alt.ids.to.remove);
}	#	fi

  #   strip leading snails from ids.  If we do not this will cause
  #       problemd below when using the S4 record ShortRead uses to 
  #       represent fastq records
#	Commenting the line out currently, since the new version of count.spikes.R
#		does not output the snail
    #ids.to.remove <- sub("@", "", ids.to.remove);

	#	This conditional is another part of error-checking, for the case when
	#		there are no reads to be removed (e.g. when both .reads.to.remove
        #		files are empty)
if(length(ids.to.remove) > 0)	{
    
    id.filter <- srFilter(function(x)   {
        !(x@id %in% ids.to.remove)
    },   #   function(x),
    name="id.filter");
    
    id.keep <- srFilter(function(x)   {
        (x@id %in% ids.to.remove)
    },   #   function(x),
    name="id.keep");

    output.fastq.records <- fastq.records[id.filter(fastq.records)];
    output.spike.records <- fastq.records[id.keep(fastq.records)]
    
} else {	#	none to remove, so the input fastq is the output fastq
    output.fastq.records <- fastq.records;
    output.spike.records <- NULL
}	#	else

num.records.removed <- length(fastq.records) - length(output.fastq.records);
cat("Removing ", num.records.removed, " fastq records\n", sep="");

     
###   write the fastq out
### Have to check for pre-existing files first
## if (file.exists(output.fastq)){
##     file.remove(output.fastq)
## }

## if (file.exists(output.spikes)){
##     file.remove(output.spikes)
## }


writeFastq(output.fastq.records, file=output.fastq, compress=FALSE, mode = "a");

# Can't write NULL fastq file, so improvise
if(is.null(output.spike.records)){
                                        #    write.table(matrix(nrow=1,ncol=1,data="No.spikes"),sep = '\t', quote = F, row.names=F, col.names=F)
    writeFastq(output.spike.records[1], file=output.spikes, compress=F, mode = "a")
} else {
    writeFastq(output.spike.records, file=output.spikes, compress=FALSE, mode = "a")
}

