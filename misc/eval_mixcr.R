###
### Evaluate MiXCR
###

library(data.table)

### Script to check MiXCR export clones results for contamination and output questionable reads to file.

### Arguments
arguments <- commandArgs(trailingOnly=T)

clone_dir_v <- arguments[1]
align_dir_v <- arguments[2]
metadata_v <- arguments[3]
fastq_info_v <- arguments[4]
qc_out_dir_v <- arguments[5]
out_dir_v <- arguments[6]

treatments_v <- c("none", "P14_spleen", "OT1_spleen", "P14+OT1_spleen")
#treatments_v <- "P14+OT1_spleen"

### Wrangling
metadata_dt <- fread(metadata_v)
fastq_info_dt <- fread(fastq_info_v); fastq_info_dt$sample <- sapply(unlist(fastq_info_dt[,1]), function(x) gsub("^.*S|\\.assembled.*", "", x))
clone_files_v <- list.files(clone_dir_v)
align_files_v <- list.files(align_dir_v)
file_name_lsv <- strsplit(clone_files_v[1], split = "S[0-9]+")
align_file_name_lsv <- strsplit(align_files_v[1], split = "S[0-9]+")

### Evaluation

qc_out_mat <- NULL

for (i in 1:length(treatments_v)){
    print(treatments_v[i])
    ## Subset tissue
    curr_meta_dt <- metadata_dt[metadata_dt$tissue == treatments_v[i],]
    for (j in 1:curr_meta_dt[,.N]){
        ## Get file
        curr_sample_v <- curr_meta_dt$sample[j]; print(curr_sample_v)
        curr_spike_type_v <- curr_meta_dt$spike_type[j]
        curr_spike_amount_v <- curr_meta_dt$spike_amount[j]
        curr_data_dt <- fread(paste0(clone_dir_v, file_name_lsv[[1]][1], "S", curr_sample_v, file_name_lsv[[1]][2]))
        curr_align_dt <- fread(paste0(align_dir_v, align_file_name_lsv[[1]][1], "S", curr_sample_v, align_file_name_lsv[[1]][2]))
        
        ## Unique Clones
        curr_num_clones_v <- curr_data_dt[,.N]
        
        ## Total Reads
        curr_reads_v <- unlist(sapply(curr_data_dt$Reads, function(x) unlist(strsplit(as.character(x), split = ',')), USE.NAMES = F))
        curr_mixcr_reads_v <- length(curr_reads_v)
        curr_total_reads_v <- fastq_info_dt[fastq_info_dt$sample == curr_sample_v,3]
        curr_complete_reads_v <- fastq_info_dt[fastq_info_dt$sample == curr_sample_v,2]
        curr_mixcr_pct_v <- round(curr_mixcr_reads_v / curr_total_reads_v * 100, digits = 2)

        ## V and J distribution
        curr_v_v <- gsub("\\*00", "", curr_data_dt$`Best V hit`)
        curr_j_v <- gsub("\\*00", "", curr_data_dt$`Best J hit`)
        curr_vj_v <- paste(curr_v_v, curr_j_v, sep = "_")

        ## Extract Read IDs
        curr_read_ids_v <- curr_align_dt[curr_align_dt$`Read id` %in% curr_reads_v, `Description R1`]

        ## QC Output
        curr_row_v <- c(curr_sample_v, treatments_v[i], curr_spike_type_v, curr_spike_amount_v, curr_num_clones_v, curr_complete_reads_v,
                        curr_total_reads_v, curr_mixcr_reads_v, curr_mixcr_pct_v,
                        length(unique(curr_v_v)), length(unique(curr_j_v)), length(unique(curr_vj_v)))
        qc_out_mat <- rbind(qc_out_mat, curr_row_v)

        ## Reads Output
        write.table(curr_read_ids_v, paste0(out_dir_v, file_name_lsv[[1]][1], "S", curr_sample_v, "_read_ids.txt"),
                    sep = '\t', quote = F, row.names = F)
    } # for j
} # for i

colnames(qc_out_mat) <- c("Sample", "Tissue", "Spike.Type", "Spike.Amount", "Unique.Clones", "Total.Reads", "Total.Despiked.Reads", "Mixcr.Reads",
                          "Mixcr.Pct", "V.segments", "J.segments", "VJ.combos")

write.table(qc_out_mat, paste0(qc_out_dir_v, file_name_lsv[[1]][1], "more_checks.txt"),
            sep = '\t', quote = F, row.names = F)
