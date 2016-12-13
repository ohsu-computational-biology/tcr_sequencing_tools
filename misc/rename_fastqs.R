# Rename files from core. Extract treatment type and place in a vector to be used in vdjtools metadata file.

# Get command line arguments
arguments <- commandArgs(trailingOnly=TRUE);

# Directory of files to be renamed. Should either be /DNAXXXLC/fastqs_from_core/fastqs/
#                                                      or
#                                                     /DNAXXXLC/peared_fastqs/assembled/
fastq.dir <- arguments[1]

# Create list of files from directory
files.in.dir <- list.files(fastq.dir)

batch <- unlist(strsplit(files.in.dir[1], "_"))[1]
  
# Create vector to contain treatment strings
#treatment <- numeric(length(files.in.dir))
sample <- numeric(length(files.in.dir))
  
  for (i in 1:length(files.in.dir)){
    curr.file <- files.in.dir[i]
    splitter <- unlist(strsplit(curr.file, "_"))
    if (identical(splitter[2], "Control")){
      treatment[i] <- paste(splitter[2], splitter[3])
      sample[i] <- splitter[4]
      new.file <- paste(splitter[1], splitter[4], splitter[5], splitter[6], sep = '_')
      file.rename(from=paste(fastq.dir, curr.file, sep='/'), paste(fastq.dir, '/', new.file, sep=''))
    } # fi
    else {
#    treatment[i] <- splitter[2]     # Takes treatment string and adds to a vector for later use
    sample[i] <- splitter[2]
    new.file <- paste(splitter[1], splitter[3], splitter[4], splitter[5], sep = '_')
    file.rename(from=paste(fastq.dir, curr.file, sep = '/'), paste(fastq.dir, '/', new.file, sep=''))
    } # else
  } # for
  
#output <- data.frame("Treatment" = treatment, "Sample" = sample)

output <- data.frame("Sample" = sample)

write.table(output,
            file=paste(batch, "treatment_codes.txt", sep = "_"),
            quote=FALSE,
            sep='\t',
            row.names=FALSE)
