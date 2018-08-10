#!/usr/bin/Rscript

extractRecord <- function(string_v, record_v, names_v){
  # Extract the number and percent values from mixcr records
  # string_v - character vector to distinguish which element of record to look at
  # record_v - assemblage of character records read in from mixcr record
  # names_v - character vectors to use as column names
  grep_v <- grep(string_v, record_v, value = T)
  if (length(grep_v) == 0){
    final_v <- c(NA, NA)
    names(final_v) <- names_v
  } else {
    partial_v <- unlist(strsplit(trimws(strsplit(grep_v, ":")[[1]][2]), " "))
    final_v <- as.numeric(gsub("\\(|%|\\)", "", partial_v))
    names(final_v) <- names_v
  }
  return(final_v)
}

