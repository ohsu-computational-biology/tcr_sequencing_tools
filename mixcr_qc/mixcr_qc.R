# This script takes an alignment and an assembly export file from MiXCR and extracts Read IDs
# of reads that successfully aligned, as well as creating a QC output.

#################
### Arguments ###
#################

arguments <- commandArgs(trailingOnly = T)

align.file <- arguments[1]
assemble.file <- arguments[2]
align.output <- arguments[3]
assemble.output <- arguments[4]
qc.output <- arguments[5]

##############
### Set Up ###
##############

# Read in Data
align.data <- read.table(align.file, sep = "\t", header = T, na.strings = c('', ' '))
assemble.data <- read.table(assemble.file, sep = "\t", header = T, na.strings = c('', ' '))

# Extract Sample Name
sample.number <- gsub(".*_S|_alignment.*", '', align.file)

# Extract Batch Name
batch <- gsub(".*DNA|LC.*", '', align.file)

# Store total aligned reads and total assembled reads as variables
total.aligned <- length(unique(align.data$Description.R1))
total.assembled <- length(align.data[complete.cases(align.data$Clone.Id),1])

#####################################
### Read IDs of Aligned Sequences ###
#####################################

# Extract read IDs for all aligned reads
aligned.ids <- as.data.frame(align.data$Description.R1)

# Create output file name
align.output.name <- paste("S", sample.number, "_align_read_ids.txt", sep = '')

# Write to an output file
write.table(aligned.ids, file = file.path(align.output, align.output.name),
              sep = '\t',
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)

#######################################
### Read IDs of Assembled Sequences ###
#######################################

### In the align export table, the sequences that eventually get assembled into clones
### have the corresponding clone ID listed. Sequences that were not assembled have NAs here.

# Extract all sequences with Clone IDs from align file
assembled.seqs <- align.data[(complete.cases(align.data$Clone.Id)),]

# Extract ids of these sequences into a data frame
assembled.ids <- as.data.frame(assembled.seqs$Description.R1)

# Create output file name
assemble.output.name <- paste("S", sample.number, "_assembled_read_ids.txt", sep = '')

# Write to an output file
write.table(assembled.ids, file = file.path(assemble.output, assemble.output.name),
              sep = '\t',
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)

####################################################
### Extract Read IDs of aligned without D region ###
####################################################

# How many reads are considered aligned even though they don't have D regions?
align.missing.d <- align.data[!(complete.cases(align.data$All.D.hits)),]

missing.d.ids <- as.data.frame(align.missing.d$Description.R1)

missing.d.output.name <- paste("S", sample.number, "_align_no_d_ids.txt", sep = '')

write.table(missing.d.ids, file = file.path(align.output, missing.d.output.name),
            sep = '\t',
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)



##############################
### Extract QC Information ###
##############################

# Percent of reads considered aligned even though they don't have D regions
align.percent.missing.d <- length(align.missing.d[,1]) / total.aligned * 100


# How many aligned don't have a Clone ID? Meaning they didn't get assembled
align.no.clone.id <- align.data[!(complete.cases(align.data$Clone.Id)),]
align.pct.no.clone <- length(align.no.clone.id[,1]) / total.aligned * 100


# How many aligned, have a D region, but didn't assemble?
align.have.d <- align.data[(complete.cases(align.data$All.D.hits)),]
align.have.d.not.assembled <- align.have.d[!(complete.cases(align.have.d$Clone.Id)),]
# As a percentage of all alignments with D regions
pct.have.d.not.assembled <- length(align.have.d.not.assembled[,1]) / 
  length(align.have.d[,1]) * 100
# As a percentage of all alignments - don't think this tells us anything
#pct.have.d.not.assembled.full <- length(align.have.d.not.assembled[,1]) /
#  total.aligned * 100

# How many assembled reads don't have D regions?
assemble.missing.d <- assemble.data[!(complete.cases(assemble.data$Best.D.hit)),]
assemble.percent.missing.d <- length(assemble.missing.d[,1]) / total.assembled * 100

# Construct QC table
qc.summary <- data.frame(sample.id = character(),
                         batch.id = character(),
                         aligned.reads = integer(),
                         assembled.reads = integer(),
                         pct.aligned.missing.d = integer(),
                         pct.aligned.not.assembled = integer(),
                         pct.aligned.have.d.not.assembled.of.aligned.w.d = integer(),
                         #pct.aligned.have.d.not.assembled.of.total.aligned = integer(),
                         pct.assembled.missing.d = integer(),
                         stringsAsFactors=FALSE)

qc.summary[1,]$sample.id <- paste("S", sample.number, sep = '')
qc.summary[1,]$batch.id <- batch
qc.summary[1,]$aligned.reads <- total.aligned
qc.summary[1,]$assembled.reads <- total.assembled
qc.summary[1,]$pct.aligned.missing.d <- align.percent.missing.d
qc.summary[1,]$pct.aligned.not.assembled <- align.pct.no.clone
qc.summary[1,]$pct.aligned.have.d.not.assembled.of.aligned.w.d <- pct.have.d.not.assembled
qc.summary[1,]$pct.assembled.missing.d <- assemble.percent.missing.d

# Create file name
qc.name <- paste("S", sample.number, "_mixcr_qc.txt", sep = '')

# Write table
write.table(qc.summary, file = file.path(qc.output, qc.name),
            quote = FALSE,
            sep = '\t',
            row.names = FALSE)

