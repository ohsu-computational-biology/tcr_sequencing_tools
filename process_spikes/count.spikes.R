##  This script:
#   

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

    #   create outputs
    output.table <- spike.table;
    output.table$spike.count <- 0;
    records.to.remove.ids <- character();

    #   count spikes - note there are two cases
    for(i in 1:length(spikes))  {
        if((i %% 10) == 0)  {
            cat("Processing spike ", i, " out of ", length(spikes), "\n", sep="");
        }   #   fi
        current.spike <- spikes[i];
        #   count spike occurences
        current.spike.counts <- vcountPattern(current.spike, sread(fastq.reads));
        output.table[i,]$spike.count <- sum(as.logical(current.spike.counts));
        cat("\t", output.table[i,]$spike.count, " spikes detected\n", sep="");
        #   record read ids (as they'll need removed later)
        records.to.remove <- which(current.spike.counts > 0);
        if(length(records.to.remove) > 0)   {
            temp.reads <- fastq.reads[records.to.remove];
            records.to.remove.ids <- c(records.to.remove.ids, 
                                        as.character(temp.reads@id));
        }   #   fi
        if(spike.length == 9)
            break;
    }   #   for i
    #   construct summary for QC purposes
    qc.summary <- data.frame(sample.id=character(),
                             num.reads=integer(),
                             num.spiked.reads=integer(),
                             pct.spiked.reads=numeric(),
                             stringsAsFactors=FALSE);
    qc.summary[1,]$sample.id <- input.fastq;
    qc.summary[1,]$num.reads <- num.fastqs;
    qc.summary[1,]$num.spiked.reads <- sum(output.table$spike.count);
    qc.summary[1,]$pct.spiked.reads <- (qc.summary$num.spiked.reads / num.fastqs) * 100;
    #   modify output.table so it'll become easily incorporated into qc.summary
    spike.count <- output.table$spike.count;
    qc.spike <- data.frame(spike.count);
    rownames(qc.spike) <- output.table$SPIKE_ID;
    qc.spike <- t(qc.spike);
    rownames(qc.spike) <- NULL;
    qc.summary <- cbind(qc.summary, qc.spike); 
    
    #   build names of output files
    #   strip ".fastq" from file name.  Not necessary, just more hygienic 
    filename.no.extension <- sub("[.][^.]*$", "", basename(input.fastq));
    count.table <- paste(output.dir, filename.no.extension, "spike.counts.", spike.length, "bp.txt", sep="");
	#	TODO:  strip .fastq away before appending ".reads.to.remove.txt"
    reads.to.remove.list <- paste(output.dir, filename.no.extension, ".reads.to.remove.txt", sep="");
    qc.file <- paste(output.dir, filename.no.extension, ".qc.txt", sep="");
    cat("Writing spike count file to: ", count.table, "\n", sep="");
    cat("Writing list of reads to remove to: ", reads.to.remove.list, "\n", sep="");
    cat("Writing QC summary to: ", qc.file, "\n", sep="");
    #   write outputs
    write.table(output.table,
                file=count.table,
                quote=FALSE,
                sep=",");

    write.table(records.to.remove.ids,
                file=reads.to.remove.list,
                quote=FALSE,
                row.names=FALSE,
                col.names=FALSE);

    write.table(qc.summary,
                file=qc.file,
                quote=FALSE,
                sep=",",
                row.names=FALSE);

}   #   count.spikes()

