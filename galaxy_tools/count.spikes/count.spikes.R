#' Count Spikes
#'
#' This function takes a fastq file and counts ocurrences of spike sequences from spike file. Generates three types of output:
#' \enumerate{
#' \item Quality control summary
#' \item A list of fastq read ids that contain spikes
#' \item A table containing spikes and counts for each spike.
#' }
#' @param input.fastq Single-end or merged Paired-end fastq file.
#' @param spike.file Tab-separated file containing spike sequences.
#' @param spike.length Any integer. Generally 9 for spike removal and 25 for normalization.
#' @param output.dir Directory in which output files will be written

#   load depdencies
suppressMessages(source("https://bioconductor.org/biocLite.R", echo=FALSE, verbose=FALSE))
suppressMessages(library(ShortRead));
suppressMessages(library(data.table))


input.fastq <- commandArgs(trailingOnly=T)[1]
spike.file <- commandArgs(trailingOnly=T)[2]
spike.length <- commandArgs(trailingOnly=T)[3]
output.count <- commandArgs(trailingOnly=T)[4]
output.remove <- commandArgs(trailingOnly=T)[5]
output.qc <- commandArgs(trailingOnly=T)[6]

###   Read in fastq file
fastq.reads <- readFastq(input.fastq);
num.fastqs <- length(fastq.reads);
cat(num.fastqs, " fastq reads to process\n", sep="");

###   open spike file
spike.table <- fread(spike.file)

###   extract spike sequences into a vector.
spikes <- spike.table$SPIKE;

###   Verify that all spikes are same length
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
    ##   count spike occurences
    ##   searches each fastq read for current spike, returns TRUE (found) or FALSE (not found)
    current.spike.counts <- vcountPattern(current.spike,
                                          sread(fastq.reads),
                                          max.mismatch = 1,
                                          with.indels = TRUE);
    
    ##   Count TRUE occurrences, add to output table, and update log
    output.table[i,]$spike.count <- sum(as.logical(current.spike.counts));
    ##   update progress
    cat("\t", output.table[i,]$spike.count, " spikes detected\n", sep="");
    
    ##   record readids (as they'll need removed later; we use the readid to do so)
    records.to.remove <- which(current.spike.counts > 0);
    if(length(records.to.remove) > 0)   {
        temp.reads <- fastq.reads[records.to.remove];
        records.to.remove.ids <- c(records.to.remove.ids,
                                   as.character(temp.reads@id));
      }   #   fi
          
    ##   In the case where spike.length == 9, all spikes are identical (since the
    ##       first nine characters of each spike are identical.  We use this break
    ##       to short-circuit the "for" loop, since we get all the information we
    ##       need in the first pass
    ##   If the first 9 characters are not identical, stop.
    if(spike.length == 9){
        check <- spikes[1] == spikes
        check <- which(check == FALSE)
        if (length(check) == 0){
            break;
        } else {
            stop("First nine characters of spikes are not identical.")
        }   #   else
    }   #   fi
            
    ##   update progress to the user
    if((i %% 10) == 0)  {
        cat("Processing spike ", i, " out of ", length(spikes), "\n", sep="");
    }   #   fi
}   #   for i

### Construct summary for QC purposes
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
    

###   write outputs
write.table(output.table,
            file=output.count,
            quote=FALSE,
            sep="\t");
  
write.table(records.to.remove.ids,
            file=output.remove,
            quote=FALSE,
            row.names=FALSE);

write.table(qc.summary,
            file=output.qc,
            quote=FALSE,
            sep='\t',
            row.names=F)

