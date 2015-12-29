
#	Assumptions:
#
#		1.  Naming format for raw clone count files:
#
#				SAMPLEX_exported_clones.txt
#
#		2.  Naming format for processed clone count files 
#			(clones with frameshifts and stop codons removed):
#
#				SAMPLEX_exported_clones_no_fssc.txt
#
#	Import required libraries
#
library(stringr);

arguments <- commandArgs(trailingOnly=TRUE);
path.to.raw.counts <- arguments[1];
path.to.processed.counts <- arguments[2];


	#	Get lists of files
	raw.counts <- list.files(path.to.raw.counts);
	processed.counts <- list.files(path.to.processed.counts);

	#	Check for same number of files
	if(length(raw.counts) != length(processed.counts))	{
		stop("Mismatch between number of raw clone counts and processed clone counts");
	}	#	fi

	#	check for parallelism of samples
	sample.id.raw.counts <- character(length(raw.counts));
	sample.id.processed.counts <- character(length(processed.counts));
	for(i in 1:length(raw.counts))	{
		sample.id.raw.counts[i] <- str_split(raw.counts[i], "_")[[1]][1];
		sample.id.processed.counts[i] <- str_split(processed.counts[i], "_")[[1]][1];
	}	#	for i

	#	Check for parallelism of files
	sample.comparison <- sample.id.raw.counts == sample.id.processed.counts;
	sample.comparison <- which(sample.comparison == FALSE);
	if(length(sample.comparison) > 0)	{
		stop("Mismatch between raw.counts and processed.counts\n");
	}	#	fi

	samples <- character();

	#	Calculate difference in file sizes
	lc.raw.counts <- NULL;
	lc.processed.counts <- NULL;
	for(i in 1:length(raw.counts))	{
		curr.raw.system.call <- paste("wc -l ", 
									  path.to.raw.counts, raw.counts[i], 
									  " | awk '{print $1}'",
									  sep="");
		curr.processed.system.call <- paste("wc -l ", 
									  path.to.processed.counts, processed.counts[i],
									  " | awk '{print $1}'",
									  sep="");
		lc.raw.counts[i] <- as.numeric(system(curr.raw.system.call, intern=TRUE));
		lc.processed.counts[i] <- as.numeric(system(curr.processed.system.call, intern=TRUE));

		samples[i] <- basename(processed.counts[i]);

		#	report status, periodically
		if((i %% 10) == 0)	{
			cat("Processing file ", i, " of ", length(raw.counts), "\n", sep="");
			cat("Current raw clone counts:  ", raw.counts[i], "\n", sep="");
			cat("Current processed clone counts: ", processed.counts[i], "\n\n", sep="");
		}	#	fi 
	}	#	for i

	#	Compare results
	records.removed <- (lc.raw.counts - lc.processed.counts);
    percent.original.clones.retained <- round(((lc.processed.counts / lc.raw.counts) * 100), digits=0);
	result.df <- data.frame(samples, lc.raw.counts, lc.processed.counts, records.removed, percent.original.clones.retained);

    #   sort data frame, to make results more intuitive
    result.df <- result.df[order(result.df$percent.original.clones.retained, decreasing=TRUE),];
	#	rename columns
	names(result.df)[1] <- "Sample";
	names(result.df)[2] <- "Clone count (original)";
	names(result.df)[3] <- "Clone count (after removal of frameshifts and stop codons)";
	names(result.df)[4] <- "Clones removed";
	names(result.df)[5] <- "Percent of original clones retained";

    #   write out results
	output.file <- "frameshift.stop.remover.QC.results.txt";
	write.table(result.df, 
				file=output.file,
				quote=FALSE,
				sep=",",
				row.names=FALSE);










