#  This script counts occurrences of spikes in a fastq file.
#
#   The script generates three types of output:
#
#       1.  Quality control summary
#       2.  A list of fastq read ids; viz., reads containing spikes
#       3.  A table containing spikes and counts for each spike
#
#   Note that runtimes are substantially shorter when spike.length
#       equals 9
#
#   load depdencies
.libPaths("/mnt/lustre1/CompBio/lib/R/library")
library(ShortRead);
library(stringr);
library(Biostrings);

arguments <- commandArgs(trailingOnly=TRUE);
input.fastq <- arguments[1];
spike.file <- arguments[2]; # for example: text_barcodesvj.txt
spike.length <- arguments[3]; # typically 9 for spike removal and 25 for normalization
output.dir <- arguments[4]; # results are output here


# Begin count.spikes function here

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


    #    trim to length (regardless of whether we're using the spike or its reverse-compl.
    spikes <- strtrim(spikes, spike.length);

    #   create outputs
    output.table <- spike.table;
    output.table$spike.count <- 0;
    records.to.remove.ids <- character();

    #   count spikes - note there are two cases
    for(i in 1:length(spikes))  {
       
        #   get a single spike 
        current.spike <- spikes[i];

        #   count spike occurences
        current.spike.counts <- vcountPattern(current.spike, sread(fastq.reads), max.mismatch = 1, with.indels = TRUE);
        output.table[i,]$spike.count <- sum(as.logical(current.spike.counts));
        cat("\t", output.table[i,]$spike.count, " spikes detected\n", sep="");

        #   record read ids (as they'll need removed later; we use the readid to do so)
        records.to.remove <- which(current.spike.counts > 0);
        if(length(records.to.remove) > 0)   {
            temp.reads <- fastq.reads[records.to.remove];
            records.to.remove.ids <- c(records.to.remove.ids, 
                                        as.character(temp.reads@id));
        }   #   fi
        
        #   In the case where spike.length == 9, all spikes are identical (since the
        #       first nine characters of each spike are identical.  We use this break
        #       to short-circuit the "for" loop, since we get all the information we
        #       need in the first pass
        #   TODO:  what if spikes don't share a common first nine characters?  Modify
        #       the code to deal with this case
        if(spike.length == 9)
            break;

        #   update progress to the user 
        if((i %% 10) == 0)  {
            cat("Processing spike ", i, " out of ", length(spikes), "\n", sep="");
        }   #   fi
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
    count.table <- paste(output.dir, filename.no.extension, ".spike.counts.", spike.length, "bp.txt", sep="");
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


