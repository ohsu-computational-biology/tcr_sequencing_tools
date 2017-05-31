###
### Test Fastq Files for Clones
###

### Take a fastq file as input, then sequentially search for different manifestations of the two 9-bp barcodes and divide reads accordingly

### Dependencies
library(data.table)
library(ShortRead)

### Arguments

arguments <- commandArgs(trailingOnly = T)

fastq_dir_v <- arguments[1]
spike_ref_v <- arguments[2]
metadata_v <- arguments[3]
out_dir_v <- arguments[4]
type_v <- arguments[5]
operation_v <- arguments[6] # "both" - search for all barcode forms and output not only hits, but also remainder of input file
                            # "remove" - only search for smallest barcode, but do not output. Only output remainder of input file
                            # "hits.only" - search for all barcode forms and only output hits

# fastq_dir_v <- "~/Desktop/OHSU/tcr_spike/data/LIB170111LC/assembled_reads/"
# #fastq_dir_v <- "~/Desktop/OHSU/tcr_spike/data/LIB170111LC/spikes/"
# spike_ref_v <- "~/Desktop/OHSU/tcr_spike/all_barcodes.txt"
# metadata_v <- "~/Desktop/OHSU/tcr_spike/data/LIB170111LC/metadata.txt"
# out_dir_v <- "~/Desktop/OHSU/tcr_spike/data/LIB170111LC/spike_search_results/"
# type_v <- "mixcr" # or despiked
# operation_v <- "remove"

### Inputs

fastq_files_v <- list.files(fastq_dir_v)
spike_ref_dt <- fread(spike_ref_v)
metadata_dt <- fread(metadata_v)


### Miscellaneous arguments

fwd_barcode_v <- "TAAGTCGAC"
rev_barcode_v <- "GTCGACTTA"

### Prepare barcodes
makeBarcodeDT <- function(barcode_v, n){
  out_barcode_dt <- as.data.table(cbind(c(barcode_v, sapply(1:n, function(x) paste0(unlist(strsplit(barcode_v, split = ''))[1:(nchar(barcode_v)-x)], collapse = '' ))),
                                        c(barcode_v, sapply(1:n, function(x) paste0(unlist(strsplit(barcode_v, split = ''))[(x+1):nchar(barcode_v)], collapse = '')))))
  colnames(out_barcode_dt) <- c("from.end", "from.start")
  return(out_barcode_dt)
}

fwd_barcode_dt <- makeBarcodeDT(fwd_barcode_v, 5)
rev_barcode_dt <- makeBarcodeDT(rev_barcode_v, 5)

### Combine into list
all_barcodes_dt <- cbind(fwd_barcode_dt, rev_barcode_dt)
all_barcodes_lsv <- as.list(as.data.table(t(all_barcodes_dt))); names(all_barcodes_lsv) <- as.character(1:length(all_barcodes_lsv)-1)
all_barcodes_lsv <- lapply(all_barcodes_lsv, unique)
if (operation_v %in%  "remove") all_barcodes_lsv <- all_barcodes_lsv[[length(all_barcodes_lsv)]]

### Function to Search Fastq File for Sequence
searchFastqForSeq <- function(fastq_data_SRQ, search_string_v){
  # Return short read object of all reads that contain search string.
  search_results_v <- vcountPattern(search_string_v, sread(fastq_data_SRQ))
  positive_indeces_v <- which(search_results_v > 0)
  positive_hits_SRQ <- fastq_data_SRQ[positive_indeces_v]
  #positive_ids_v <- as.character(positive_hits_SRQ@id)
  return(positive_hits_SRQ)
}

### Function to combine Fastq Files
combineFastq <- function(fastq_data_lsSRQ){
  first_fastq_SRQ <- fastq_data_lsSRQ[[1]]
  for (i in 2:length(fastq_data_lsSRQ)){
    first_fastq_SRQ <- ShortRead::append(first_fastq_SRQ, fastq_data_lsSRQ[[i]])
  }
  return(first_fastq_SRQ)
}

# search_results_lsls <- list()
# 
# for (j in 1:fwd_barcode_dt[,.N]){
#   curr_name_v <- as.character(j-1)
#   curr_barcode_v <- unique(unlist(c(fwd_barcode_dt[j,], rev_barcode_dt[j,])))
#   for (k in 1:curr_barcode_v){
#     for (i in 1:length(fastq_files_v)){
#       curr_fastq_path_v <- paste0(fastq_dir_v, fastq_files_v[i])
#       curr_search_results_v <- searchFastqForSeq(curr_fastq_path_v, curr_barcode_v[k])
#     } # for i
#   } # for k
# } # for j

### Make empty variables to hold information
barcode_summary_mat <- matrix(nrow = length(fastq_files_v), ncol = (2+fwd_barcode_dt[,.N]))
barcode_search_data_lslsv <- list()
if (operation_v %in% c("both", "remove")) no_hits_lsv <- list()

### Perform search for each fastq file
for (i in 1:length(fastq_files_v)){
  ## Update
  if (i %% 10 == 0) print(c("Currently on file ", i))
  ## Get a file
  curr_fastq_data_SRQ <- readFastq(file.path(fastq_dir_v, fastq_files_v[i]))
  curr_ids_v <- as.character(curr_fastq_data_SRQ@id)
  curr_summary_row_v <- c(fastq_files_v[i], length(curr_fastq_data_SRQ))
    curr_results_lsv <- list()
#    print(all_barcodes_lsv)
  for (j in 1:length(all_barcodes_lsv)){
    ## Get a set of barcodes
    if (operation_v == "remove") {
      curr_barcode_v <- all_barcodes_lsv
    } else {
      curr_barcode_v <- all_barcodes_lsv[[j]]
    }
#      print(operation_v)
#      print(all_barcodes_lsv)
#      print(curr_barcode_v)
    ## Return a list of reads with each barcode, one element for each barcode, element contains short read object of corresponding reads
      curr_search_results_lsSRQ <- sapply(curr_barcode_v, function(x) searchFastqForSeq(curr_fastq_data_SRQ, x))
#      print(c("Currently on file ", i))
#      print("Test")
#      print(str(curr_search_results_lsSRQ))
#      print(class(curr_search_results_lsSRQ))
#      print(c("Length of search results ", length(curr_search_results_lsSRQ)))
    curr_search_results_SRQ <- combineFastq(curr_search_results_lsSRQ)
    ## Extract IDs out in a similarly-structued list
    # old way
    #curr_search_results_ids_lsv <- lapply(curr_search_results_lsSRQ, function(x) as.character(x@id))
    #curr_unique_ids_v <- unlist(curr_search_results_ids_lsv); curr_unique_ids_v <- unique(curr_unique_ids_v)
    curr_search_results_ids_v <- as.character(curr_search_results_SRQ@id)
    ## Get all unique IDs
    curr_unique_ids_v <- unique(curr_search_results_ids_v)
    ## Get indeces of all unique reads and first occurrence of duplicate reads
    curr_indeces_v <- match(curr_unique_ids_v, curr_search_results_ids_v)
    ## Subset short read object to only contain these sequences
    curr_search_results_SRQ <- curr_search_results_SRQ[curr_indeces_v]
    ## Find indeces in original fastq of positive hits, and remove them from original for next iteration!
    curr_keep_ids_v <- curr_ids_v[!(curr_ids_v %in% curr_unique_ids_v)]
    curr_keep_indeces_v <- match(curr_keep_ids_v, curr_ids_v)
    curr_fastq_data_SRQ <- curr_fastq_data_SRQ[curr_keep_indeces_v]
    curr_ids_v <- as.character(curr_fastq_data_SRQ@id)
    ## Update summary row and output
    curr_summary_row_v <- c(curr_summary_row_v, length(curr_unique_ids_v))
    if (operation_v == "remove"){
      curr_results_lsv[[ "5" ]] <- curr_search_results_SRQ
    } else {
      curr_results_lsv[[ names(all_barcodes_lsv)[j] ]] <- curr_search_results_SRQ
    } # fi
    
  } # for j
  if (operation_v %in% c("both", "remove")) no_hits_lsv[[ fastq_files_v[i] ]] <- curr_fastq_data_SRQ
  if (operation_v %in% c("both", "hits.only")) {
    barcode_search_data_lslsv[[ fastq_files_v[i] ]] <- curr_results_lsv
    barcode_summary_mat[i,] <- curr_summary_row_v
  }
} # for i

if (operation_v %in% c("both", "hits.only")){
  barcode_summary_dt <- as.data.table(barcode_summary_mat); colnames(barcode_summary_dt) <- c("Sample", "Total.Reads", as.character(c(0, 1, 2, 3, 4, 5)))
  change_cols_v <- as.character(c(0,1,2,3,4,5))
  barcode_summary_dt[,(change_cols_v) := lapply(.SD, as.numeric), .SDcols = change_cols_v]
  barcode_summary_dt$sample.num <- as.integer(gsub("^.*_S|_assemb.*|\\.assemb.*", "", barcode_summary_dt$Sample))
  barcode_summary_dt <- merge(barcode_summary_dt, metadata_dt, by.x = "sample.num", by.y = "sample")
  
  write.table(barcode_summary_dt, file = paste0(out_dir_v, type_v, "_search_summary.txt"),
              sep = '\t', quote = F, row.names = F)
}


if (operation_v %in% c("both", "hits.only"))
for (i in 1:length(barcode_search_data_lslsv)){
  curr_name_v <- names(barcode_search_data_lslsv)[i]
  curr_data_lsv <- barcode_search_data_lslsv[[i]]
  for (j in 1:length(curr_data_lsv)){
    curr_barcode_name_v <- names(curr_data_lsv)[j]
    if (length(curr_data_lsv[[j]]) > 0){
      writeFastq(curr_data_lsv[[j]], file = file.path(out_dir_v, curr_barcode_name_v, curr_name_v), compress = F)
    } # fi
  } # for j
} # for i

if (operation_v %in% c("both", "remove")){
  for (i in 1:length(no_hits_lsv)){
    curr_name_v <- names(no_hits_lsv)[i]
    curr_data_v <- no_hits_lsv[[i]]
    if (length(curr_data_v) > 0){
      writeFastq(curr_data_v, file = file.path(out_dir_v, "remainder", curr_name_v), compress = F)
    } # fi
  } # for i
} # fi
