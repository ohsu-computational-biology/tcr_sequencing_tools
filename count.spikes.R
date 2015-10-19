##  This script removes spiked records from a fastq file

#   load depdencies
library(ShortRead);
library(stringr);
library(Biostrings);

count.spikes <- function(input.fastq, 
                        spike.file,
                        spike.length=34, 
                        output.dir,
                        direction="FWD")  {

    #   Read in fastq file
    fastq.reads <- readFastq(input.fastq);
    num.fastqs <- length(fastq.reads);
    cat(num.fastqs, " fastq reads to process\n", sep="");

    #   open and prepare spike file
    spike.table <- read.csv(spike.file,
                            stringsAsFactors=FALSE,
                            sep=" ");
    spikes <- spike.table$SPIKE;
    #   TODO:  add some error-checking, verifying that all spikes are of same length, etc.

    #   modify spikes as required
    if(direction == "REV")  {
        #   first convert to a DNAString, to convert to reverse-complement,
        #       then back to a character
        #   I do not believe that DNAString is vectorized, so we loop
       for(i in 1:length(spikes))   {
           spikes[i] <- as.character(reverseComplement(DNAString(spikes[i]))); 
       }   #   for i 
    }   #   fi
    #    trim to length (regardless of whether we're using the spike or its reverse-compl.
    spikes <- strtrim(spikes, spike.length);

    #   create output
    output.table <- spike.table;
    output.table$spike.count <- 0;

    #   count spikes - note there are two cases
    for(i in 1:length(spikes))  {
        if((i %% 10) == 0)  {
            cat("Processing spike ", i, " out of ", length(spikes), "\n", sep="");
        }   #   fi
        current.spike <- spikes[i];
        #   count spike occurences
        current.spike.counts <- vcountPattern(current.spike, sread(fastq.reads));
        output.table[i,]$spike.count <- sum(as.logical(current.spike.counts));
        #   record read ids (as they'll need removed later)
    }   #   for i

    #   build sames of output files
    count.table <- paste(output.dir, basename(input.fastq), ".counts.", spike.length, "bp.txt", sep="");
    reads.to.remove.list <- paste(output.dir, basename(input.fastq), ".reads.to.remove.txt", sep="");
    cat("Writing spike count file to: ", count.table, "\n", sep="");
    cat("Writing list of reads to remove to: ", reads.to.remove.list, "\n", sep="");
    
    #   write output
    write.table(output.table,
                file=count.table,
                quote=FALSE,
                sep=",");

}   #   count.spikes()

