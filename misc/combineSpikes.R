### Aggregate spike counts together

library(data.table)

arguments <- commandArgs(trailingOnly = T)

spike_dir_v <- arguments[1]
out_dir_v <- arguments[2]

spike_files_v <- list.files(spike_dir_v)

total_data_mat <- matrix(nrow=260, ncol = length(spike_files_v))
names_v <- character(length = length(spike_files_v))

batch_v <- gsub("LC_S.*", "", spike_files_v[1])

for (i in 1:length(spike_files_v)){
    curr_file_v <- spike_files_v[i]
    curr_name_v <- strsplit(curr_file_v, split = "_|\\.")[[1]][2]
    names_v[i] <- curr_name_v
    
    curr_data_dt <- fread(paste0(spike_dir_v, curr_file_v))

    total_data_mat[,i] <- curr_data_dt$V6
}

total_data_mat <- cbind(curr_data_dt$V4, curr_data_dt$V5, total_data_mat)


colnames(total_data_mat) <- c("V", "J", names_v)

write.table(total_data_mat, file = paste0(out_dir_v, batch_v, "_25bp_counts.txt"),
            sep = '\t', quote = F, row.names = F)
