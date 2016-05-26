#   Script that assesses quality of the "remove spikes" step in the TCRseq
#       pipeline

# Get command line arguments
arguments <- commandArgs(trailingOnly=TRUE);
path.to.raw.fastqs <- arguments[1];		# /DNAXXXXLC/peared_fastqs/assembled/
path.to.despiked.fastqs <- arguments[2];	# /DNAXXXXLC/mixcr/despiked_fastqs



	#	Get lists of files
	raw.fastqs <- list.files(path.to.raw.fastqs);
	despiked.fastqs <- list.files(path.to.despiked.fastqs);

	#	Strip file names to raw sample names
	raw.fastq.samples <- sub("[.][^.]*$", "", raw.fastqs); 
	raw.fastq.samples <- sub("[.][^.]*$", "", raw.fastq.samples); 
	despiked.fastq.samples <- sub("[.][^.]*$", "", despiked.fastqs); 
	despiked.fastq.samples <- sub("[.][^.]*$", "", despiked.fastq.samples); 
	despiked.fastq.samples <- sub("[.][^.]*$", "", despiked.fastq.samples); 
	despiked.fastq.samples <- sub("[.][^.]*$", "", despiked.fastq.samples); 

	#	Check for parallelism of files
	sample.comparison <- raw.fastq.samples == despiked.fastq.samples;
	sample.comparison <- which(sample.comparison == FALSE);
	if(length(sample.comparison) > 0)	{
		stop("Mismatch between raw fastq samples and despiked fastq samples\n");
	}	#	fi

	#	Calculate difference fastq files
	lc.raw.fastqs <- NULL;
	lc.despiked.fastqs <- NULL;
	for(i in 1:length(raw.fastqs))	{
		curr.raw.system.call <- paste("wc -l ", 
									  path.to.raw.fastqs, raw.fastqs[i], 
									  " | awk '{print $1}'",
									  sep="");
		curr.despiked.system.call <- paste("wc -l ", 
									  path.to.despiked.fastqs, despiked.fastqs[i],
									  " | awk '{print $1}'",
									  sep="");
		lc.raw.fastqs[i] <- as.numeric(system(curr.raw.system.call, intern=TRUE));
		lc.despiked.fastqs[i] <- as.numeric(system(curr.despiked.system.call, intern=TRUE));

		#	report status, periodically
		if((i %% 10) == 0)	{
			cat("Processing file ", i, " of ", length(raw.fastqs), "\n", sep="");
			cat("Current raw fastq:  ", raw.fastqs[i], "\n", sep="");
			cat("Current despiked fastq: ", despiked.fastqs[i], "\n\n", sep="");
		}	#	fi 
	}	#	for i

	#	Compare results
    #   We divide by four (in a number of cases) because each fastq record contains four lines
	records.removed <- (lc.raw.fastqs - lc.despiked.fastqs);
	records.removed <- records.removed / 4;
	lc.raw.fastqs <- lc.raw.fastqs / 4;
	lc.despiked.fastqs <- lc.despiked.fastqs / 4;
    percent.original.reads.retained <- (lc.despiked.fastqs / lc.raw.fastqs);
    percent.original.reads.retained <- round(percent.original.reads.retained * 100, digits=1);
	result.df <- data.frame(despiked.fastqs, 
                            lc.raw.fastqs, 
                            lc.despiked.fastqs,
                            records.removed,
                            percent.original.reads.retained);
    #   sort data frame, to make results more intuitive
    result.df <- result.df[order(result.df$percent.original.reads.retained, decreasing=TRUE),];
    #   rename columns
    names(result.df)[1] <- "Sample (despiked)";
    names(result.df)[2] <- "Number of reads (original)";
    names(result.df)[3] <- "Number of reads (after removal of spiked reads";
    names(result.df)[4] <- "Number of reads removed";
    names(result.df)[5] <- "Percent of original reads retained";
    #   write out results
	output.file <- "remove.spikes.QC.result.txt";
    write.table(result.df,
                file=output.file,
                quote=FALSE,
                sep=",",
                row.names=FALSE);

