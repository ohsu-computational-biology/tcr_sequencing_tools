#	RScript that calculates number of unique clonotypes, Shannon entropy, and clonality for multiple clonotype files
#
#   Input is a directory containing one or more files with MiXCR output format
#   
#   Load necessary libraries
library(entropy);
library(data.table);
library(tcR);
library(ggplot2);

#	Get command-line arguments
arguments <- commandArgs(trailingOnly=TRUE);
#   Directory should contain all and only .txt files from exportClones, post-normalization with the following
#   format:
#     Clone count   Clone fraction    Clonal sequence(s)    AA. Seq. CDR3   Best V Hit    Best J Hit    V segments
#     J segments    Normalized clone count    Normalized clone fraction

cloneDir_v <- arguments[1];    # Typically .../dhaarini/DNAXXXXLC/normalization/normalized_clones/
metadata_v <- arguments[2]     # Typically $data/QC/meta
outDir_v <- arguments[3];      # Typeically $data/QC/homeo
countOutDir_v <- arguments[4]  # Typically $data/QC/cloneDiv
toWrite_v <- arguments[5]      # true or false to write out count files and pdfs or not
tissue_v <- arguments[6]       # specify a specific tissue to subset by. If NA, will use all tissues
type_v <- arguments[7]         # sometimes divide by a certain category of treatment, rather than each treatment. If NA, will not divide

### Get files and names
cloneFiles_v <- list.files(cloneDir_v)
cloneFiles_v <- cloneFiles_v[order(as.numeric(gsub("^.*_S|_.*$", "", cloneFiles_v)))]
cloneNames_v <- sapply(cloneFiles_v, function(x) strsplit(x, split = "_")[[1]][2], USE.NAMES=F)
batchName_v <- strsplit(cloneFiles_v[1], split = "_")[[1]][1]

### Read in data
cloneData_lsdt <- sapply(cloneFiles_v, function(x) fread(file.path(cloneDir_v, x)), simplify = F)
names(cloneData_lsdt) <- cloneNames_v

metadata_dt <- fread(metadata_v)

### Get appropriate columns and subset

tissueCol_v <- grep("issue", colnames(metadata_dt), value = T)
sampleCol_v <- grep("ample", colnames(metadata_dt), value = T)[1]

if (!is.na(tissue_v)) {
    ## Subset metadata
    metadata_dt <- metadata_dt[get(tissueCol_v) == tissue_v,]
    ## Subset clone data
    keepSamp_v <- paste0("S", metadata_dt[, get(sampleCol_v)])
    cloneData_lsdt <- cloneData_lsdt[keepSamp_v]
} #fi

### Get treatments to run analysis on
if (!is.na(type_v)){
    treatCol_v <- type_v
} else {
    treatCol_v <- grep("eatment", colnames(metadata_dt), value = T)
} # fi

treatments_v <- unique(metadata_dt[, get(treatCol_v)])

### Initialize empty variables
finalHomeo_mat <- NULL
finalRows_v <- NULL

if (toWrite_v){
    pdf(file = paste0(outDir_v, batchName_v, "_homeostasis.pdf"))
} #fi

for (i in 1:length(treatments_v)){

    ## Get Treatment
    currTreat_v <- treatments_v[i]

    ## Subset files
    currFiles_v <- paste0("S", metadata_dt[get(treatCol_v) == currTreat_v, get(sampleCol_v)]) # LIB170728LC
    currData_lsdt <- cloneData_lsdt[currFiles_v]

    ## Print nrow
    nrow_v <- sapply(currData_lsdt, nrow)
    print(c("Treatment" = currTreat_v, nrow_v))

    ## Remove zero-count clones
    currData_lsdt <- lapply(currData_lsdt, function(x) x[`Normalized clone fraction` > 0,])

    ## Run clonal space homeostasis
    currHomeo_mat <- clonal.space.homeostasis(.data = currData_lsdt, .prop.col = "Normalized clone fraction")
    currHomeo_dt <- data.table("Sample" = rownames(currHomeo_mat),
                               currHomeo_mat)

    ## Check that frequencies sum to 1
    cumFreq_v <- unname(apply(currHomeo_mat, 1, function(x) round(sum(x), digits = 10)))
    bad_v <- which(cumFreq_v != 1)
    if (length(bad_v) > 0) {
        stop("Incorrect cumulative frequency")
    } # fi

    ## Write counts
    if (toWrite_v){
        write.table(currHomeo_dt, file = paste0(countOutDir_v, currTreat_v, "_cloneDiv.txt"),
                    sep = '\t', quote = F, row.names = F)
    } # fi

    ## Take mean of matrix
    currMean_v <- apply(currHomeo_mat, 2, mean)

    ## Add to list
    finalRows_v <- c(finalRows_v, currTreat_v)
    finalHomeo_mat <- rbind(finalHomeo_mat, currMean_v)

    ## Make plot
    currHomeo_gg <- vis.clonal.space(currHomeo_mat) +
        ggtitle(paste0(currTreat_v, " Clonal Space Homeostasis")) +
        theme(plot.title = element_text(hjust = 0.5))

    if (toWrite_v){
        print(currHomeo_gg)
    } #fi
    
} # for i

### Finish final table construction
rownames(finalHomeo_mat) <- finalRows_v
finalHomeo_dt <- data.table("Treatment" = finalRows_v, finalHomeo_mat)

### Write
if (toWrite_v){
    ## Table
    write.table(finalHomeo_dt, file = paste0(countOutDir_v, batchName_v, "_cloneDiv.txt"),
                row.names = F, quote = F, sep = '\t')
    ## Plot
    finalHomeo_gg <- vis.clonal.space(finalHomeo_mat) +
        ggtitle(paste0("Treat-wise Clonal Space Homeostasis")) +
        theme(plot.title = element_text(hjust = 0.5))

    print(finalHomeo_gg)
    dev.off()
} #fi
