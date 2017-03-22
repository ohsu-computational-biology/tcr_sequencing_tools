# Rename files from core. 
# Extract treatment type and place in a vector to be used in vdjtools metadata file.
# Most common format: DNAXXXXLC_S[0-9]+_S[0-9]+_R[1-2]_001.fastq.
# Want to remove extraneous S[0-9]+ 
# There may be an extra treatment or tissue sample identifier. This is currently set up to work with 1 extra identifier

# Get command line arguments
arguments <- commandArgs(trailingOnly=TRUE);

# Directory of files to be renamed. Should either be /DNAXXXLC/fastqs_from_core/fastqs/
#                                                      or
#                                                     /DNAXXXLC/peared_fastqs/assembled/

fastq.dir <- arguments[1]
out.dir <- arguments[2]
extra_segment_pos <- arguments[3]

# Fix if no output directory specified
if (is.na(out.dir)){
   out.dir <- "./"
   message("No output directory given. Writing to current directory.")
} # fi

# Create list of files from directory
files.in.dir <- list.files(fastq.dir)

# Get batch identifier
batch <- unlist(strsplit(files.in.dir[1], "_"))[1]
  
# Add vector for extra segment
if (!is.na(extra_segment_pos)){
	extra_segment <- numeric(length(files.in.dir))
}

# Add vector for samples
sample <- numeric(length(files.in.dir))
  
  for (i in 1:length(files.in.dir)){

    # Get file and split
    curr.file <- files.in.dir[i]
    split.file <- unlist(strsplit(curr.file, "_"))

    # Old catch for control samples...may not work as is, so if there are "Control" samples, please check this!
    if (identical(split.file[2], "Control")){
      treatment[i] <- paste(split.file[2], split.file[3])
      sample[i] <- split.file[4]
      new.file <- paste(splitfile[1], split.file[4], split.file[5], split.file[6], sep = '_')
      file.rename(from=paste(fastq.dir, curr.file, sep='/'), paste(fastq.dir, '/', new.file, sep=''))
    } # fi
    else {

    # Save extra segment, if present
    if (!is.na(extra_segment_pos)){
	extra_segment[i] <- split.file[extra_segment_pos]
    } # fi

    # Save one of the samples for output
    sample[i] <- split.file[2]
    # combine the rest to a new file name
    new.file <- paste(split.file[1],		# batch
			 split.file[3], 	# sample number
			 split.file[4], 	# R1/R2
			 split.file[5], 	# 001.fastq
			 sep = '_')
    # Rename file
    file.rename(from=paste(fastq.dir, curr.file, sep = '/'), paste(fastq.dir, '/', new.file, sep=''))

    } # fi
  } # for

if (!is.na(extra_segment_pos)){
   output <- data.frame("Sample" = sample, "ExtraCol" = extra_segment)
} else {
   output <- data.frame("Sample" = sample)
} # fi

write.table(output,
            file=paste0(out.dir, batch, "_treatment_codes.txt"),
            quote=FALSE,
            sep='\t',
            row.names=FALSE)
