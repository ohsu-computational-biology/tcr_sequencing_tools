#!/usr/bin/env Rscript

### DIVERSITY ANALYSIS

###	RScript that calculates number of unique clonotypes, Shannon entropy, and clonality for multiple clonotype files
### Optional - run the calculations on subsets of the top clones instead of the entire repertoire

####################
### DEPENDENCIES ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

suppressMessages(library(data.table))
suppressMessages(library(tcR))
suppressMessages(library(optparse))

#################
### FUNCTIONS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################

dsName <- function(dataSetName_v){
  # Return last directory from file path and file name of galaxy dataset.dat path
  # dataSetName_v - single element of a galaxy dataset path (e.g. /database/files/000/dataset123.dat)
  splitName_v <- strsplit(dataSetName_v, split = "/")[[1]]
  x_v <- length(splitName_v)
  newName_v <- paste(splitName_v[c((x_v-1),x_v)], collapse = "/")
  return(newName_v)
}

my.entropy.plugin <- function (freqs, unit = c("log", "log2", "log10")) {
  ### Taken from 'entropy' package. Can't install package because of dependency issues in conda
  unit = match.arg(unit)
  freqs = freqs/sum(freqs)
  H = -sum(ifelse(freqs > 0, freqs * log(freqs), 0))
  if (unit == "log2") 
    H = H/log(2)
  if (unit == "log10") 
    H = H/log(10)
  return(H)
}

####################
### COMMAND LINE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

optlist <- list(
  make_option(
    c("-f", "--files"),
    type = "character",
    help = "Normalized clonotype files. If specifying the freq divisions argument (--divisions), this file must have a column
    which specifies those divisions. Standard normalized clone files don't have this, must be added by Group Clones tool."
  ),
  make_option(
    c("-n", "--names"),
    type = "character",
    help = "Names of normalized clonotype files"
  ),
  make_option(
    c("-o", "--out"),
    type = "character",
    help = "Name of output file"
  ),
  make_option(
    c("-l", "--old"),
    type = "logical",
    help = "TRUE - input has hold column names (or have un-normalized data). FALSE - use updated names."
  ),
  make_option(
    c("-d", "--divisions"),
    type = "character",
    help = "Comma-separated list of either integers or frequency divisions, with no spaces. Determines which subsets of the clonotype repertoire will be analyzed.
    If blank, the entire reperotire will be analyzed. List of integers will do 'top set', so '10,25,50' will run analysis on top10, top25, and top50 clones.
    list of frequency divisions will do 'freq group', so 'Small,Medium,Large' will analyzed clones in the small, medium, and large clonal frequency groups."
  )
)

p <- OptionParser(usage = "%prog -f files -n names -o out -l old -d divisions",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

cloneFiles_v <- args$files
cloneNames_v <- args$names
outFile_v <- args$out
old_v <- args$old
divisions_v <- args$divisions

### Split command-line arguments
cloneFiles_v <- unlist(strsplit(cloneFiles_v, ','))

### Read in names
cloneNames_dt <- fread(cloneNames_v, header = F)

### Handle divisions
print("Raw divisions")
print(divisions_v)

if (!is.null(divisions_v)){ # 2018-10-05 changed is.na to is.null because switched from commandArgs() to optparse.
  ## Split
  divisions_v <- strsplit(divisions_v, split = ',')[[1]]
  ## Check for numeric
  divCheck_v <- as.numeric(divisions_v)
  if (is.na(divCheck_v[1])) {
    print("Clone frequency divisions set.")
  } else {
    print("Top clone divisions set.")
    divisions_v <- divCheck_v
  } # fi 
} else {
  divisions_v <- "noDiv"
} # fi


# if (!is.na(divisions_v)) {
#   if (!is.na(freq_v)) {
#     stop("Both frequency group divisions and top clone divisions are set! Only one may be set. Please unset one or the other.")
#   } else {
#     divisions_v <- strsplit(divisions_v, split = ",")[[1]]
#     divisions_v <- as.numeric(divisions_v)
#     print("Top clone divisions set.")
#   } # fi
# } else {
#   if (!is.na(freq_v)){
#     divisions_v <- strsplit(freq_v, split = ",")[[1]]
#     print("Clone frequency divisions set.")
#   } else {
#     divisions_v <- "noDiv"
#     print("No clone subset.")
#   } # fi
# } # fi
print("Processed divisions")
print(divisions_v)

############
### BODY ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############

for (j in 1:length(divisions_v)) {
  
  ## Get divisions
  currDiv_v <- divisions_v[j]
  
  ## Create arrays
  cloneFileNames_v <- sampleNames_v <- character(length(cloneFiles_v))
  calculated.entropies <- unique.clones <- clonality <- gini <- true <- numeric(length(cloneFiles_v))
  max.clonal.freq <- max.clone.count <- norm.entropy <- adaptive.clonality <- clonality
  
  ## Update
  cat(sprintf("Working on division %s\n", currDiv_v))
  
  for(i in 1:length(cloneFiles_v))	{
    ##
    ## DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##
    
    ## Get a file, the sample name, and the database name
    currClone_v <- cloneFiles_v[i];
    sampleNames_v[i] <- cloneNames_dt$V1[i]
    cloneFileNames_v[i] <- dsName(currClone_v)
    
    ## Read in data
    currData_dt <- fread(currClone_v)
    
    ## Skip if empty
    if (nrow(currData_dt) == 0) next
    
    ## Get column names
    if (old_v) {
      hasNorm_v <- grep("Normalized", colnames(currData_dt), value = T)
      if (length(hasNorm_v) > 0){
        column_v <- "Normalized clone fraction"
        count_v <- "Normalized clone count"
      } else {
        column_v <- grep("Clone fraction|cloneFraction", colnames(currData_dt), value = T)
        count_v <- grep("Clone count|cloneCount", colnames(currData_dt), value = T)
      } # fi
    } else {
      column_v <- "nb.clone.fraction"
      count_v <- "nb.clone.count"
    } # fi
    
    ## If taking divisions, need to sort and subset
    if (currDiv_v != "noDiv") {
      if (is.numeric(currDiv_v)){
        currData_dt <- currData_dt[order(currData_dt[[column_v]], decreasing = T),]
        currData_dt <- currData_dt[1:currDiv_v]
      } else {
        currData_dt <- currData_dt[Div == currDiv_v,]
      } # fi
    } # fi
    
    ##
    ## CALCULATIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##
    
    ## UNIQUE CLONES
    unique.clones[i] <- currData_dt[,.N]
    
    ## ENTROPY
    #calculated.entropies[i] <- entropy(currData_dt[[column_v]], method="ML", unit="log");
    calculated.entropies[i] <- my.entropy.plugin(freqs = currData_dt[[column_v]], unit = "log")
    
    ## CLONALITY
    clonality[i] <- 1 - (calculated.entropies[i] / log(unique.clones[i]))
    
    ## NORMALIZED SHANNON ENTROPY
    norm.entropy[i] <- calculated.entropies[i] / log(unique.clones[i])
    
    ## ADAPTIVE CLONALITY (inverse of normalized entropy)
    adaptive.clonality[i] <- 1 / norm.entropy[i]
    
    ## GINI
    gini[i] <- gini(.data = currData_dt[[column_v]], .do.norm = F)
    true[i] <- diversity(.data = currData_dt[[column_v]], .do.norm = T)
    
    ## Change clone freq column to a percentage
    currData_dt[[column_v]] <- currData_dt[[column_v]] * 100
    
    ## MAX CLONAL FREQ
    max.clonal.freq[i] <- round(max(currData_dt[[column_v]]), digits = 4)
    
    ## MAX CLONE COUNT
    max.clone.count[i] <- max(currData_dt[[count_v]])
  } # for i
  
  ###   create output data.frame
  output.df <- data.frame(cloneFileNames_v, sampleNames_v, calculated.entropies, norm.entropy, unique.clones, clonality,
                          adaptive.clonality, max.clonal.freq, max.clone.count, gini, true);
  
  colnames(output.df) <- c("File", "Sample", "Shannon Entropy", "Normalized Entropy", "Unique Clonotypes", "Clonality",
                           "Adaptive Clonality", "Max Clonal Freq", "Max Clone Count", "Gini_Index", "True_Diversity")
  
  ### If running divisions, add the division to the column names
  if (currDiv_v != "noDiv") colnames(output.df) <- paste(colnames(output.df), currDiv_v, sep = "_")
  
  ### Combine into final data frame
  if (j == 1) {
    final_df <- output.df
  } else {
    final_df <- cbind(final_df, output.df)
  }
  
} # for j

###
### If ran divisions, will have duplicates of the file/sample columns
###

### This will extract all other columns (but removes all copies of file and sample)
keepCols_v <- grep("File_|Sample_", colnames(final_df), value = T, invert = T)

### This will add back the 1st versions of File and sample
addCols_v <- grep("_noDiv", paste(c("File", "Sample"), divisions_v[1], sep = "_"), value = T, invert = T)

### Extract these columns
final_df <- final_df[,c(addCols_v, keepCols_v)]

###   write output
write.table(final_df, 
            file=outFile_v,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
