####################################################
### Remove monoclonal contamination from samples ###
####################################################

library(data.table)

# We have three monoclonal sequences that have been used at various points throughout the project. After their introduction, 
# we have noticed an overabundance of their presence in later samples. This script searches through clone files exported by
# Mixcr exportClones using the V identity, J identity, and AA sequence of the three clonal contaminants and records their 
# fastq read id for later removal.

### Clones

# p14
  # V133, J2-4, CASSDAGGRNTLYF
# ot1
  # V121, J2-7, CASSRANYEQYF
# el4
  # V15, J2-3, CASSTGTETLYF

##############
### Set Up ###
##############

arguments <- commandArgs(trailingOnly = T)

clone.dir <- arguments[1]        # #data/mixcr/export_clones
align.dir <- arguments[2]
out.dir <- arguments[3]          # $data/qc/decontam_id

### Sort files
clone.files <- list.files(clone.dir)
clone.files <- clone.files[order(as.numeric(gsub(".*_S|_align.*|_clones.txt", '', clone.files)))]

align.files <- list.files(align.dir)
align.files <- align.files[order(as.numeric(gsub(".*_S|_align.*|_clones.txt", '', align.files)))]

### Batch name for output files
batch_v <- unlist(strsplit(clone.files[1], split = "_"))[1]

### Empty variable to hold read ids
outputHeaders_v <- NULL

###################
### Calculation ###
###################


# Iterate through each clone file, removing monoclonal contaminants
for (i in 1:length(clone.files)){
    
    ## Sample name
    temp <- unlist(strsplit(clone.files[i], split = "_"))
    currName_v <- grep("S[0-9]+", temp, value = T)
    tempCheck <- unlist(strsplit(align.files[i], split = "_"))
    checkName_v <- grep("S[0-9]+", tempCheck, value = T)

    ## Stop if names don't match
    if (checkName_v != currName_v) stop("Incorrect matching of clone and align files")

    ## Update
    print(currName_v)

    ## Sample number
    currNum_v <- gsub("S", "", currName_v)
    
    ## Get files
    curr.clone <- fread(paste(clone.dir, clone.files[i], sep = ''))
    curr.align <- fread(paste(align.dir, align.files[i], sep = ''), select = c("readId", "descrR1"))

    ## Get column names
    idCol_v <- grep("Clone ID|cloneId", colnames(curr.clone), value = T)[1]
    rawCountCol_v <- grep("Clone count|cloneCount", colnames(curr.clone), value = T)
    vCol_v <- grep("Best V hit|bestVHit", colnames(curr.clone), value = T)
    jCol_v <- grep("Best J hit|bestJHit", colnames(curr.clone), value = T)
    readCol_v <- grep("Reads|reads", colnames(curr.clone), value = T)
    seqCol_v <- grep("aa.*CDR3|AA.*CDR3", colnames(curr.clone), value = T)
    fracCol_v <- grep("Clone fraction|cloneFraction", colnames(curr.clone), value = T)
    print(c(idCol_v, rawCountCol_v, vCol_v, jCol_v, readCol_v, seqCol_v, fracCol_v))

    ## Add index for ranking
    curr.clone$index <- seq(1, length(curr.clone[[idCol_v]]))

    ## Update column names and contents for V and J segments
    curr.clone$`V segments` <- gsub("TRB|\\*00", "", curr.clone[[vCol_v]])
    curr.clone$`V segments` <- gsub("-", "", curr.clone$`V segments`)
    curr.clone$`J segments` <- gsub("TRB|\\*00", "", curr.clone[[jCol_v]])

    
    ## Extract clones (take just the first entry if more than one clone matches the criteria)

    ##
    ## p14
    ##
    offending.clone.1 <- curr.clone[(curr.clone$`V segments` == "V133" # p14
        & curr.clone$`J segments` == "J2-4" &
          curr.clone[[seqCol_v]] == "CASSDAGGRNTLYF"),]
    
    offending.clone.1 <- offending.clone.1[1,] # subset for first only

    ## Get read IDs
    offending.clone.1.ids <- unlist(strsplit(offending.clone.1[[readCol_v]], split = ","))

    ## Get fastq ids
    offending.clone.1.headers <- curr.align[readId %in% offending.clone.1.ids, descrR1]

    ## Add to global list
    outputHeaders_v <- c(outputHeaders_v, offending.clone.1.headers)

    ##
    ## ot1
    ##
    offending.clone.2 <- curr.clone[(curr.clone$`V segments` == "V121" # ot1
        & curr.clone$`J segments` == "J2-7" &
          curr.clone[[seqCol_v]] == "CASSRANYEQYF"),]

    offending.clone.2 <- offending.clone.2[1,] # subset for first only

    ## Get read IDs
    offending.clone.2.ids <- unlist(strsplit(offending.clone.2[[readCol_v]], split = ","))

    ## Get fastq ids
    offending.clone.2.headers <- curr.align[readId %in% offending.clone.2.ids, descrR1]

    ## Add to global list
    outputHeaders_v <- c(outputHeaders_v, offending.clone.2.headers)

    ##
    ## el4
    ##
    offending.clone.3 <- curr.clone[(curr.clone$`V segments` == "V15" # el4
        & curr.clone$`J segments` == "J2-3" &
          curr.clone[[seqCol_v]] == "CASSTGTETLYF"),]
    offending.clone.3 <- offending.clone.3[1,] # subset for first only
  
    ## Get read IDs
    offending.clone.3.ids <- unlist(strsplit(offending.clone.3[[readCol_v]], split = ","))

    ## Get fastq ids
    offending.clone.3.headers <- curr.align[readId %in% offending.clone.3.ids, descrR1]

    ## Add to global list
    outputHeaders_v <- c(outputHeaders_v, offending.clone.3.headers)

    ##
    ## output
    ##

    ## Format output
    outputHeaders_dt <- as.data.table(outputHeaders_v)
    colnames(outputHeaders_dt) <- "Reads"

    ## Create output name
    out.name <- file.path(out.dir, paste(batch_v, currName_v, "contamIDs.txt", sep = "_")) 
  
    write.table(outputHeaders_dt, out.name,
                    sep = '\t', quote = F, row.names = F)

    ## Update progress
    if (i %% 10 == 0){print(c("Currently on clone ", i))}
    
} # for

