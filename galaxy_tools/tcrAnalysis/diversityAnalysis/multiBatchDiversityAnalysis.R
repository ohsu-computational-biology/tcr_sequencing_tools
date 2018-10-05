### Analysis Script

### Calculate a variety of summary and diversity statistics for all of the normalized clonotype files

### Calculate Shannon Entropy, normalized Shannon Entropy, clonality, adaptive's clonality, unique clones, max clone count, gini index, true diversity

####################
### DEPENDENCIES ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

suppressMessages(library(data.table))
suppressMessages(library(tcR))

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
    c("-b", "--batch"),
    type = "integer",
    default = 1,
    help = "Index of batch ID from sample name, when split by '_'."
  )
)

p <- OptionParser(usage = "%prog -f files -n names -o out -l old -d divisions",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

cloneFiles_v <- args$files
cloneNames_v <- args$names
batchIndex_v <- args$batch
outFile_v <- args$out

### Get all of the files and name them
cloneFiles_v <- unlist(strsplit(cloneFiles_v, ','))
cloneNames_v <- fread(cloneNames_v, header = F)$V1
names(cloneFiles_v) <- cloneNames_v
print(cloneFiles_v)

### Get batches and order them
batches_v <- unique(sapply(cloneNames_v, function(x) strsplit(x, split = "_")[[1]][batchIndex_v]))
batches_v <- sort(batches_v)

### Sort files by batch and then by sample within batch
orderedNames_v <- unname(unlist(sapply(batches_v, function(x) {
	y <- cloneNames_v[grep(x, cloneNames_v)]
	z <- y[order(as.numeric(gsub(".*_S|_.*$|\\..*$}", '', y)))]
	return(z)}, simplify = F)))
cloneFiles_v <- cloneFiles_v[orderedNames_v]

############
### BODY ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############

### Create empty arrays
fileNames_v <- galaxyNames_v <- batchCol_v <- character(length(cloneFiles_v))
calculated.entropies <- unique.clones <- clonality <- gini <- true <- numeric(length(cloneFiles_v))
max.clonal.freq <- max.clone.count <- norm.entropy <- adaptive.clonality <- clonality

for(i in 1:length(cloneFiles_v))	{

    ###
    ### DATA
    ###
  cat(sprintf("on iteration: %s\n", i))

    ## Get a clone file to process and read it in
    currFile_v <- cloneFiles_v[i];
    currData_dt <- fread(currFile_v)
    
    ## Get the galaxy name and the actual name
    galaxyNames_v[i] <- currGalaxy_v <- dsName(currFile_v)
    fileNames_v[i] <- currName_v <- names(cloneFiles_v)[i]
    
    ## Get batch
    currBatch_v <- strsplit(currName_v, split = "_")[[1]][batchIndex_v]
    batchCol_v[i] <- currBatch_v

    ## Get columns
    column_v <- "normFreq"
    count_v <- "normC"

    ###
    ### CALCULATIONS
    ###

    ## Count number of lines in file, i.e. number of unique clonotypes
    unique.clones[i] <- currData_dt[,.N]

    ## Calculate entropy
    calculated.entropies[i] <- my.entropy.plugin(currData_dt[[column_v]], unit = "log")
    
    ## Calculate clonality
    clonality[i] <- 1 - (calculated.entropies[i] / log(unique.clones[i]))

    ## Calculate clonality as the inverse of normalized shannon entropy
    ## Normalized shannon entropy
    norm.entropy[i] <- calculated.entropies[i] / log(unique.clones[i])
    ## New clonality
    adaptive.clonality[i] <- 1 / norm.entropy[i]

    ## tcR Gini and True Diversity
    gini[i] <- gini(.data = currData_dt[[column_v]], .do.norm = F)
    true[i] <- diversity(.data = currData_dt[[column_v]], .do.norm = F)
	
    ## Change clone frequency column to a percentage
    currData_dt[[column_v]] <- currData_dt[[column_v]] * 100
    
    ## Calculate Max. clonotype frequency
    max.clonal.freq[i] <- round(max(currData_dt[[column_v]]), digits = 4)

    ## Record maximum clone count
    max.clone.count[i] <- max(currData_dt[[count_v]])

    ##   update progress
    if((i %%10) == 0)   {
        cat("Processing file ", i, " (out of ", length(cloneFiles_v), ")\n", sep="");
    } # fi

} # for i

### Create output data.frame
output.df <- data.frame(galaxyNames_v, batchCol_v, fileNames_v, calculated.entropies, norm.entropy, unique.clones, clonality,
                        adaptive.clonality, max.clonal.freq, max.clone.count, gini, true)

colnames(output.df) <- c("File", "Batch", "Sample", "Shannon Entropy", "Normalized Entropy", "Unique Clonotypes", "Clonality",
                         "Adaptive Clonality", "Max Clonal Freq", "Max Clone Count", "Gini_Index", "True_Diversity")

print(head(output.df))

### Write output
write.table(output.df, 
            file=outFile_v,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
