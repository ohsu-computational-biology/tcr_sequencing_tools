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

    #   modify spikes as required
    

    #   create output
    output.table <- spike.table;
    output.table$spike.count <- 0;

    #   count spikes - note there are two cases
#    for(i in 1:length(spikes))  {
    for(i in 1:10)  {
        if((i %% 10) == 0)  {
            cat("Processing spike ", i, " out of ", length(spikes), "\n", sep="");
        }   #   fi
        if(direction == "FWD")  {   #   use the spike as-is
            #   TODO:  add some error-checking, e.g. with file name
            current.spike <- spikes[i];
        }   else {  #   use the reverse-complement of the spike in the spike file 
            current.spike <- as.character(reverseComplement(DNAString(spikes[i])));
            #   TODO:  add some error-checking, e.g. with file name
        }   #   else

        #   count spike occurences
        current.spike.counts <- vcountPattern(current.spike, sread(fastq.reads));
        output.table[i,]$spike.count <- sum(as.logical(current.spike.counts));

        #   record read ids (as they'll need removed later)

        
    }   #   for i

    #   build sames of output files
    count.table <- paste(output.dir, basename(spike.file), ".counts.txt", sep="");
    reads.to.remove.list <- paste(output.dir, basename(spike.file), ".reads.to.remove.txt", sep="");
    cat("Writing spike count file to: ", count.table, "\n", sep="");
    cat("Writing list of reads to remove to: ", reads.to.remove.list, "\n", sep="");
    
    #   write output
    write.table(output.table,
                file=count.table,
                quote=FALSE,
                sep=",");

}   #   count.spikes()

