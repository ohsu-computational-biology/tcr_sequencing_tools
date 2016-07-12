# Both options produce exactly the same output.

# Input Files Option 1
#arguments <- commandArgs(trailingOnly=TRUE);
#goo <- as.numeric(arguments[1])
#coo <- as.numeric(arguments[2])

# Output Files Option 1
#out.file <- arguments[3]
#out.file2 <- arguments[4]

# Input Files Option 2
#goo <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
#coo <- as.numeric(commandArgs(trailingOnly = TRUE)[2])

# Output Files Option 2
#out.file <- commandArgs(trailingOnly = TRUE)[3]
out.file2 <- commandArgs(trailingOnly = TRUE)[1]

# Calculations with input
#tot <- goo * coo
#foo <- matrix(seq(1:tot), nrow=goo, ncol=coo)
boo <- matrix(seq(1:10), nrow=1, ncol=10)
colnames(boo) <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")

# Write outputs
#write.table(foo, sep='\t', file = out.file, row.names=F)
write.table(boo, sep='\t', file = out.file2, col.names=T, row.names=F)