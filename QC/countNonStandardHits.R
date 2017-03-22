###
### Count Non-Standard Hits
###

### MiXCR can align to regions that we do not specifically amplify for. We would like to summarize the amounts of different regions.

### Dependencies
require(data.table)

### Arguments
arguments <- commandArgs(trailingOnly = T)
clone_dir_v <- arguments[1] #  data/normalization/decontam OR data/normalization/normalized_clones
spike_file_v <- arguments[2]
out_dir_v <- arguments[3]

### Set Up
clone_files_v <- list.files(clone_dir_v)
spike_data_dt <- fread(spike_file_v)
standard_j_regions_v <- unique(spike_data_dt$J)
standard_v_regions_v <- unique(spike_data_dt$V)
standard_v_regions_v <- gsub("-", '', standard_v_regions_v)
standard_v_regions_v <- c(standard_v_regions_v, "V121", "V122")
standard_v_regions_v <- standard_v_regions_v[-which(standard_v_regions == "V1212")]

nonstandard_v_lsdt <- list()
nonstandard_j_lsdt <- list()

for (i in 1:length(clone_files_v)){
    ## Get data
    curr_file_v <- clone_files_v[i]
    curr_name_v <- strsplit(curr_file_v, split = "_")[[1]][2]
    curr_data_dt <- fread(paste0(clone_dir_v, curr_file_v))

    ## Subset for non-standard only
    curr_nonstandard_v_v <- curr_data_dt[!(curr_data_dt$`V segments` %in% standard_v_regions_v), `V segments`]
    curr_nonstandard_j_v <- curr_data_dt[!(curr_data_dt$`J segments` %in% standard_j_regions_v), `J segments`]

    ## Count Occurrences
    if (length(curr_nonstandard_v_v) > 0){
        curr_nonstandard_v_dt <- as.data.table(table(curr_nonstandard_v_v))
    } else {
        curr_nonstandard_v_dt <- NULL
    }

    if (length(curr_nonstandard_j_v) > 0){
        curr_nonstandard_j_dt <- as.data.table(table(curr_nonstandard_j_v))
    } else {
        curr_nonstandard_j_dt <- NULL
    }

    ## Add to summary
    nonstandard_v_lsdt[[curr_name_v]] <- curr_nonstandard_v_dt
    nonstandard_j_lsdt[[curr_name_v]] <- curr_nonstandard_j_dt
} # for i

## Get all non-standard J regions
all_nonstandard_v_v <- unique(unlist(lapply(nonstandard_v_lsdt, function(x) x$curr_nonstandard_v_v)))
all_nonstandard_j_v <- unique(unlist(lapply(nonstandard_j_lsdt, function(x) x$curr_nonstandard_j_v)))

## Create output
if (!is.null(all_nonstandard_v_v)){
    nonstandard_v_output_mat <- matrix(nrow = length(nonstandard_v_lsdt), ncol = length(all_nonstandard_v_v))
    colnames(nonstandard_v_output_mat) <- all_nonstandard_v_v
    rownames(nonstandard_v_output_mat) <- names(nonstandard_v_lsdt)
    for (i in 1:length(nonstandard_v_lsdt)){
        curr_dt <- nonstandard_v_lsdt[[i]]
        curr_name_v <- names(nonstandard_v_lsdt)[i]
        for (j in curr_dt[,.N]){
            curr_v_v <- curr_dt$curr_nonstandard_v_v[j]
            curr_count_v <- curr_dt$N[j]
            nonstandard_v_output_mat[rownames(nonstandard_v_output_mat) == curr_name_v,
                                     colnames(nonstandard_v_output_mat) == curr_v_v] <- curr_count_v
        } # for j
    } # for i
    nonstandard_v_output_dt <- as.data.table(nonstandard_v_output_mat)
    nonstandard_v_output_dt$sample <- rownames(nonstandard_v_output_mat)
    nonstandard_v_output_dt <- nonstandard_v_output_dt[,c(ncol(nonstandard_v_output_dt), 1:(ncol(nonstandard_v_output_dt)-1)), with = F]
    write.table(nonstandard_v_output_mat, file = paste0(out_dir_v, "nonstandardVHits.txt"), quote = F, sep = '\t', row.names = F)
} else {
    print("There were no nonstandard V regions in this batch!")
} #fi
