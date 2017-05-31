###
### Aggregate clone files from a batch into a single file. If possible, include metadata information as well.
###

### Argument Descriptions
     ## clone_dir_v - direcotry of decontaminated, normalized clone files
     ## out_dir_v - directory to write final output file to
     ## metadata_v - optional. columns of treatment and other information to be added to files. 
          ## If any files in clone_dir_v are not mentioned in metadata, they will be removed. 
          ## Must have these columns:
               ## sample - numeric sample number
               ## files - full file name as it appears in clone_dir_v
               ## any other information, such as treatment, input concentration, etc.

### Dependencies
library(data.table)

### Arguments
arguments <- commandArgs(trailingOnly = T)

clone_dir_v <- arguments[1]
out_dir_v <- arguments[2]
metadata_v <- arguments[3]


#metadata_v <- "~/Desktop/OHSU/tcr_spike/data/LIB170111LC/temp"
#metadata_v <- "~/Desktop/OHSU/tcr_spike/data/blood_analysis/170303_meta.txt"
#metadata_v <- NULL
#clone_dir_v <- "~/Desktop/OHSU/tcr_spike/data/LIB170111LC/export_clones/"
#clone_dir_v <- "~/Desktop/OHSU/tcr_spike/data/LIB170303LC/normalized_clones/"
#out_dir_v <- "~/Desktop/OHSU/tcr_spike/data/LIB170111LC/test_aggregate/"
#out_dir_v <- "~/Desktop/OHSU/tcr_spike/data/LIB170303LC/aggregate/test/"

### Get clone files, sort by number, and extract base file name components (i.e. prefix and suffix).
clone_files_v <- list.files(clone_dir_v)
clone_files_v <- clone_files_v[order(as.numeric(gsub("^.*_S|_.*$|\\..*$", '', clone_files_v)))]
file_name_v <- unlist(strsplit(clone_files_v[1], split = "S[0-9]*")) # use regex in case S1 was removed for some reason upstream

### If exists, read in metadata
if (file.exists(metadata_v)){

  ## Get data
  metadata_dt <- fread(metadata_v)

  ## Remove specified files
  orig_file_num_v <- length(clone_files_v); new_file_num_v <- metadata_dt[,.N]
  print(sprintf(c("Removing %d files from consideration because not listed in metadata file.",
          "Please check that this is correct."), orig_file_num_v-new_file_num_v))
  clone_files_v <- clone_files_v[clone_files_v %in% metadata_dt$files]

  ## Check order
  metadata_dt <- metadata_dt[match(clone_files_v, metadata_dt$files),]
  if (!identical(metadata_dt$files, clone_files_v)){
    stop("Mismatch between clone files and metadata information")
  } # fi
} else {
  print("No metadata file specified. This information is important for downstream analysis based on treatment groups and other variables. Please consider specifying a file!")
}

### Create empty matrix to hold aggregated data
complete_clone_mat <- NULL

#for (i in 1:metadata_dt[,.N]){
for (i in 1:length(clone_files_v)){
  ## Update periodically
  if (i %% 10 == 0) print(sprintf("Currently on clone %d", i))

  ## Get sample number
  curr_sample_v <- gsub("^.*_S|_.*$|\\..*$", '', clone_files_v[i])

  ## Read in data, order by Normalized clone fraction, and create new id
  curr_data_dt <- fread(paste0(clone_dir_v, clone_files_v[i]))
  curr_data_dt <- curr_data_dt[order(-`Normalized clone fraction`)]
  curr_data_dt$new_id <- seq(1:curr_data_dt[,.N])

  ## Add metadata information, if applicaple
  if (file.exists(metadata_v)){
    for (j in 1:ncol(metadata_dt)){
      ## Get current value and its name
      curr_metadata_value_v <- unlist(metadata_dt[i,j, with = F])
      curr_metadata_name_v <- names(curr_metadata_value_v)
      ## Repeat it for every row in current clone file
      curr_data_dt[[curr_metadata_name_v]] <- rep(curr_metadata_value_v, curr_data_dt[,.N])
      } # for j
  } else {
    curr_data_dt$sample <- rep(curr_sample_v, curr_data_dt[,.N])
  }# fi

  ## Add to final output matrix
  complete_clone_mat <- rbind(complete_clone_mat, curr_data_dt)

} # for i

write.table(complete_clone_mat, file = file.path(out_dir_v, paste0(file_name_v[1], "aggregate_clone_info.txt")),
            sep = '\t', quote = F, row.names = F)
