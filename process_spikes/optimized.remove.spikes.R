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
library(ShortRead);
library(stringr);

#   TODO: handle warnings more better.  At least read.delim() is generating warnings
options(warn=-1);   

#remove.fastqs <- function(input.fastq, forward.reads.to.remove, reverse.reads.to.remove) {
   

	arguments <- commandArgs(trailingOnly=TRUE);
	input.fastq <- arguments[1];
	forward.reads.to.remove <- arguments[2];
	reverse.reads.to.remove <- arguments[3];

    #   Read in fastq file
    fastq.records <- readFastq(input.fastq);
    num.fastqs <- length(fastq.records);
    cat(num.fastqs, " fastq reads to process\n", sep="");

	#	create variables to store read ids to remove
	forward.ids.to.remove <- character();
	reverse.ids.to.remove <- character();
    
	#	Some samples returned an empty list of IDs to remove, raising an 
	#		exception in the original filtering code.  We added some checks here
	#		to address this possibility.  TODO:  refine this
	file.size.forward <- file.size(forward.reads.to.remove);
	if(file.size.forward > 0)	{
		#   Read in lists of fastqs to be removed and consolidate into one variable
		#   TODO:  find a cleaner way to import
		forward.ids.to.remove <- read.delim(forward.reads.to.remove,
									header=FALSE,
									stringsAsFactors=FALSE);
		forward.ids.to.remove <- forward.ids.to.remove[,1];
		alt.forward.ids.to.remove <- str_replace(forward.ids.to.remove, " 1", " 2");
		forward.ids.to.remove <- c(forward.ids.to.remove, alt.forward.ids.to.remove);
	}	#	fi

	file.size.reverse <- file.size(reverse.reads.to.remove);
	if(file.size.reverse > 0)	{
		reverse.ids.to.remove <- read.delim(reverse.reads.to.remove,
									header=FALSE,
									stringsAsFactors=FALSE);
		reverse.ids.to.remove <- reverse.ids.to.remove[,1];
		alt.reverse.ids.to.remove <- str_replace(reverse.ids.to.remove, " 2", " 1");
		reverse.ids.to.remove <- c(reverse.ids.to.remove, alt.reverse.ids.to.remove);
	}
    
    ids.to.remove <- union(forward.ids.to.remove, reverse.ids.to.remove);
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

		output.fastq.records <- fastq.records[id.filter(fastq.records)];
	}	else	{	#	none to remove, so the input fastq is the output fastq
		output.fastq.records <- fastq.records;
	}	#	else
    
	num.records.removed <- length(fastq.records) - length(output.fastq.records);
    cat("Removing ", num.records.removed, " fastq records\n", sep="");
     
    #   write the fastq out
    output.file.name <- paste(input.fastq, ".removed.fastq", sep="");
	cat("Writing output to: ", output.file.name, "\n", sep="");
    writeFastq(output.fastq.records, output.file.name, compress=FALSE);

#}   #   remove.fastqs()

