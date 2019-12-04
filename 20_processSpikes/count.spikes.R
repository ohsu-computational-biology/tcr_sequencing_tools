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

### Loading ShortRead should load BiocParallel, but I was getting a weird error and this fixed it.
library(BiocParallel, lib.loc="/home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/from_tempwork/myRPkg")
suppressMessages(suppressWarnings(library(ShortRead)));
#library(stringr);
#library(Biostrings);
suppressMessages(library(data.table));

arguments <- commandArgs(trailingOnly=TRUE);
input.fastq <- arguments[1];
spike.file <- arguments[2]; # for example: text_barcodesvj.txt
spike.length <- arguments[3]; # typically 9 for spike removal and 25 for normalization
output.dir <- arguments[4]; # results are output here


# Begin count.spikes function here

###   Read in fastq file
fastq.reads <- readFastq(input.fastq);
num.fastqs <- length(fastq.reads);
cat(num.fastqs, " fastq reads to process\n", sep="");


###   open and prepare spike file
spike.table <- fread(spike.file)

spikes <- spike.table$SPIKE;

spike_lengths <- sapply(spikes, function(x) length(x), USE.NAMES = F)

if (length(unique(spike_lengths)) != 1) warning("Not all spikes appear to be same length")


###    trim to length (regardless of whether we're using the spike or its reverse-compl.
spikes <- strtrim(spikes, spike.length);


###   create outputs
output.table <- spike.table;
output.table$spike.count <- 0;
records.to.remove.ids <- character();


###   count spikes - note there are two cases

for(i in 1:length(spikes))  {
    
    ##   get a single spike
    current.spike <- spikes[i];

    if (spike.length == 9) {
        ## Search reverse barcode as well
        reverse.spike <- "GTCGACTTA"
    } # fi
    

    ##   count spike occurences
    ##   searches each fastq read for current spike, returns TRUE (found) or FALSE (not found)
    current.spike.counts <- vcountPattern(current.spike,
                                          sread(fastq.reads),
                                          max.mismatch = 1,
                                          with.indels = TRUE);

    ## count revere spike occurrences
    if (spike.length == 9){
        reverse.spike.counts <- vcountPattern(reverse.spike,
                                              sread(fastq.reads),
                                              max.mismatch = 1,
                                              with.indels = TRUE)
        current.spike.counts <- current.spike.counts + reverse.spike.counts
    } # fi
    

    ## Count TRUE occurrences, add to output table, and update log
    output.table[i,]$spike.count <- sum(as.logical(current.spike.counts)); 
    cat("\t", output.table[i,]$spike.count, " spikes detected\n", sep="");

    ##   record read ids (as they'll need removed later; we use the readid to do so)
    ## Grab the indexes of all reads with spikes
    records.to.remove <- which(current.spike.counts > 0);

    ## If there are some, get their IDs from the fastq object
    if(length(records.to.remove) > 0)   {
        temp.reads <- fastq.reads[records.to.remove];
        records.to.remove.ids <- c(records.to.remove.ids,
                                   as.character(temp.reads@id));
        }   #   fi
        
    ##   In the case where spike.length == 9, all spikes are identical (since the
    ##   first nine characters of each spike are identical.  We use this break to
    ##   short-circuit the "for" loop, since we get all the information we
    ##   need in the first pass
    ##   If spike.length != 9, user is updated and loop repeats.
    ##   TODO:  what if spikes don't share a common first nine characters?  Modify the code to deal with this case
    if(spike.length == 9)
        break;
   
    ##   update progress to the user
    if((i %% 10) == 0)  {
        cat("Processing spike ", i, " out of ", length(spikes), "\n", sep="");
    }   #   fi
}   #   for i


###   construct summary for QC purposes
qc.summary <- data.frame(sample.id=character(),
                         num.reads=integer(),
                         num.spiked.reads=integer(),
                         pct.spiked.reads=numeric(),
                         stringsAsFactors=FALSE);
qc.summary[1,]$sample.id <- input.fastq;
qc.summary[1,]$num.reads <- num.fastqs;
qc.summary[1,]$num.spiked.reads <- sum(output.table$spike.count);
qc.summary[1,]$pct.spiked.reads <- (qc.summary$num.spiked.reads / num.fastqs) * 100;

###   modify output.table so it'll become easily incorporated into qc.summary
spike.count <- output.table$spike.count;
qc.spike <- data.frame(spike.count);
rownames(qc.spike) <- output.table$SPIKE_ID;
qc.spike <- t(qc.spike);
rownames(qc.spike) <- NULL;
qc.summary <- cbind(qc.summary, qc.spike);

### Change read vector to data.table
records.to.remove.ids.dt <- as.data.table(records.to.remove.ids)
colnames(records.to.remove.ids.dt) <- "Reads"
    
###   build names of output files
###   strip ".fastq" from file name.  Not necessary, just more hygienic
filename.no.extension <- sub("[.][^.]*$", "", basename(input.fastq));
count.table <- paste(output.dir, "counts/", filename.no.extension, ".spike.counts.", spike.length, "bp.txt", sep="");
reads.to.remove.list <- paste(output.dir, "reads_to_remove/", filename.no.extension, ".reads.to.remove.txt", sep="");
qc.file <- paste(output.dir, "qc/", filename.no.extension, ".qc.txt", sep="");

### Update
cat("Writing spike count file to: ", count.table, "\n", sep="");
cat("Writing list of reads to remove to: ", reads.to.remove.list, "\n", sep="");
cat("Writing QC summary to: ", qc.file, "\n", sep="");

###   write outputs

write.table(output.table,
            file=count.table,
            quote=FALSE,
            row.names=FALSE,
            sep="\t");


write.table(records.to.remove.ids.dt,
            file=reads.to.remove.list,
            quote=FALSE,
            row.names=FALSE);

write.table(qc.summary,
            file=qc.file,
            quote=FALSE,
            sep="\t",
            row.names=FALSE);


