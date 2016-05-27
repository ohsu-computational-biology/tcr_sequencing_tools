# Aggregate mixcr_qc files into one large data frame

#################
### Arguments ###
#################

arguments <- commandArgs(trailingOnly = T)

qc.dir <- arguments[1]


# List files in directory
qc.files <- list.files(qc.dir)

# Create empty data.frame
qc.df <- data.frame()

# Append files to it
for (i in 1:length(qc.files)){
  curr.qc.data <- read.table(file.path(qc.dir, qc.files[i]), sep = '\t', header = T,
                             stringsAsFactors = FALSE)
  qc.df <- rbind(qc.df, curr.qc.data)
} # for

batch <- curr.qc.data$batch.id

output.name <- paste(batch, "_mixcr_qc_summary.txt", sep = '')

write.table(qc.df, file=file.path(qc.dir, output.name), 
            quote = F,
            row.names = F,
            sep = '\t')
