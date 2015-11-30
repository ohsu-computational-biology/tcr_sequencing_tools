
#	Assumptions:
#
#		1.  Naming format for raw fastq files:
#
#				SAMPLE_X.assembled.fastq
#
#		2.  Naming format for despiked fastq files:
#
#				SAMPLE_X.assembled.fastq.removed.fastq

evaluate.work <- function(path.to.raw.fastqs, path.to.despiked.fastqs)	{

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
	records.removed <- ((lc.raw.fastqs - lc.despiked.fastqs) / 4);
	samples <- 1:length(raw.fastqs);
	result.df <- data.frame(samples, records.removed);
	output.file <- "results.txt";
	cat("Records removed, per sample:\n", file=output.file);
	#write(result.df, file=output.file );
	cat("\nSummary statistics on records removed:\n", file=output.file);
	print(summary(records.removed));

}	#	evaluate.work()










