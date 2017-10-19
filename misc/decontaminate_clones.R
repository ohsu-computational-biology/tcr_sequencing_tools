####################################################
### Remove monoclonal contamination from samples ###
####################################################

library(data.table)

# We have three monoclonal sequences that have been used at various points throughout the project. After their introduction, 
# we have noticed an overabundance of their presence in later samples. This script searches through clone files exported by
# Mixcr exportClones using the V identity, J identity, and AA sequence of the three clonal contaminants and removes them from
# the file. The output is the exact same as our usual MiXCR export, minus these clones

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
out.dir <- arguments[2]          # $data/normalization/decontam
qc.dir <- arguments[3]           # $data/QC
which.outputs <- arguments[4]    # "clones", "qc", or "both"
metadata.file <- arguments[5]    # path to file or "NULL" character string
keep.mono <- arguments[6]
keep.mono <- F
mono.dir <- arguments[7]

### Sort files
clone.files <- list.files(clone.dir)
clone.files <- clone.files[order(as.numeric(gsub(".*_S|_align.*", '', clone.files)))]

### Empty matrix for output summary
contam_reads_mat <- matrix(nrow = length(clone.files), ncol = 2)

### Batch name for output files
batch_v <- unlist(strsplit(clone.files[1], split = "_"))[1]

### Empty QC Matrix
contamination.qc <- matrix(nrow = length(clone.files), ncol = 14)

### Read in metadata
if (metadata.file == "NULL") {
    metadata_dt <- NA
} else {
    metadata_dt <- fread(metadata.file)
}

###################
### Calculation ###
###################


# Iterate through each clone file, removing monoclonal contaminants
for (i in 1:length(clone.files)){
    ## Counter for how many reads are devoted to contaminated sequences
    curr_contam_count_v <- 0
    
    ## Sample name
    temp <- unlist(strsplit(clone.files[i], split = "_"))
    currName_v <- grep("S[0-9]+", temp, value = T)
    #currName_v <- unlist(strsplit(clone.files[i], split = "_"))[2]

    ## Update
    print(currName_v)

    ## Sample number
    currNum_v <- gsub("S", "", currName_v)
    
    ## Get a clone file
    curr.clone <- fread(paste(clone.dir, clone.files[i], sep = ''))

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

    ## Get first QC info
    orig.unique.count <- length(curr.clone[[idCol_v]])
    orig.total.count <- sum(curr.clone[[rawCountCol_v]])
    

    ## Update column names and contents for V and J segments
    curr.clone$`V segments` <- gsub("TRB|\\*00", "", curr.clone[[vCol_v]])
    curr.clone$`V segments` <- gsub("-", "", curr.clone$`V segments`)
    curr.clone$`J segments` <- gsub("TRB|\\*00", "", curr.clone[[jCol_v]])

    
    ## Extract clones (take just the first entry if more than one clone matches the criteria)
    ## p14

    if (is.na(metadata_dt) || !("p14" %in% unlist(strsplit(metadata_dt[metadata_dt$sample == currNum_v, `mono`], split = ',')))) {
        offending.clone.1 <- curr.clone[(curr.clone$`V segments` == "V133" # p14
            & curr.clone$`J segments` == "J2-4" &
              curr.clone[[seqCol_v]] == "CASSDAGGRNTLYF"),]
    
        offending.clone.1 <- offending.clone.1[1,] # subset for first only

        ## Update count of contaminated sequences
        if (is.na(offending.clone.1[[readCol_v]][1])) {
            count.clone.1 <- 0
        } else {
            count.clone.1 <- offending.clone.1[[rawCountCol_v]]
            read.count.clone.1 <- length(unlist(strsplit(as.character(offending.clone.1[[readCol_v]]), split = ',')))
        }
        curr_contam_count_v <- curr_contam_count_v + count.clone.1
    

        ## QC
        p14.rank <- offending.clone.1$index[1]
        p14.count <- count.clone.1
    } else {
        cat("p14 not removed from sample", currName_v, "\n")
        p14.rank <- NA
        p14.count <- NA
        count.clone.1 <- NA
        read.count.clone.1 <- NA
    }
    

    ## ot1
    if (is.na(metadata_dt) || !("ot1" %in% unlist(strsplit(metadata_dt[metadata_dt$sample == currNum_v, `mono`], split = ',')))) {
        offending.clone.2 <- curr.clone[(curr.clone$`V segments` == "V121" # ot1
            & curr.clone$`J segments` == "J2-7" &
              curr.clone[[seqCol_v]] == "CASSRANYEQYF"),]
        offending.clone.2 <- offending.clone.2[1,] # subset for first only


        ## Update count of contaminated sequences
        if (is.na(offending.clone.2[[readCol_v]][1])) {
            count.clone.2 <- 0
        } else {
            count.clone.2 <- offending.clone.2[[rawCountCol_v]]
            read.count.clone.2 <- length(unlist(strsplit(as.character(offending.clone.2[[readCol_v]]), split = ',')))
        }
        curr_contam_count_v <- curr_contam_count_v + count.clone.2
    
        ## QC
        ot1.rank <- offending.clone.2$index[1]
        ot1.count <- count.clone.2
    } else {
        cat("ot1 not removed from sample", currName_v, "\n")
        ot1.rank <- NA
        ot1.count <- NA
        count.clone.2 <- NA
        read.count.clone.2 <- NA
    }


    ## el4
    if (is.na(metadata_dt) || !("el4" %in% unlist(strsplit(metadata_dt[metadata_dt$sample == currNum_v, `mono`], split = ',')))) {
        offending.clone.3 <- curr.clone[(curr.clone$`V segments` == "V15" # el4
            & curr.clone$`J segments` == "J2-3" &
              curr.clone[[seqCol_v]] == "CASSTGTETLYF"),]
        offending.clone.3 <- offending.clone.3[1,] # subset for first only
  
        if (is.na(offending.clone.3[[readCol_v]][1])) {
            count.clone.3 <- 0
            read.count.clone.3 <- 0
        } else {
            count.clone.3 <- offending.clone.3[[rawCountCol_v]]
            read.count.clone.3 <- length(unlist(strsplit(as.character(offending.clone.3[[readCol_v]]), split = ',')))
        }
        curr_contam_count_v <- curr_contam_count_v + count.clone.3
    
        ## QC
        el4.rank <- offending.clone.3$index[1]
        el4.count <- count.clone.3
    } else {
        cat("el4 not removed from sample", currName_v, "\n")
        el4.rank <- NA
        el4.count <- NA
        count.clone.3 <- NA
        read.count.clone.3 <- NA
    }
    

    ## Add count of contaminated reads to output matrix
    contam_reads_mat[i,] <- c(currName_v, curr_contam_count_v)

    
    ## Combine contaminants together into a data frame with the original clone data frame
    combo <- rbind(curr.clone, offending.clone.1, offending.clone.2, offending.clone.3)
    ## Remove duplicated entries (i.e. offending clones) from combined dataframe
    new.clone <- combo[!duplicated(combo, fromLast = T) & seq(nrow(combo)) <= nrow(curr.clone),]
    new.clone <- new.clone[!(is.na(new.clone[[idCol_v]])),]

    ## If desired, save monoclonal sequences
    if (keep.mono) {
        mono.clone <- rbind(offending.clone.1, offending.clone.2, offending.clone.3)
        write.table(mono.clone, file = paste0(mono.dir, batch_v, "_", currName_v, "_monoclonal_counts.txt"),
                    sep = '\t', quote = F, row.names = F)
    }

    ## QC Data
    remaining.count <- orig.total.count - p14.count - ot1.count - el4.count
    qc.row <- c("Sample" = currName_v, "Orig.Unique.Clones" = orig.unique.count, "Orig.total.Clones" = orig.total.count,
                "Remaining.Clones" = remaining.count, "Contam.clones" = curr_contam_count_v,
                "p14.rank" = p14.rank, "p14.count" = p14.count, "p14.reads" = read.count.clone.1,
                "ot1.rank" = ot1.rank, "ot1.count" = ot1.count, "ot1.reads" = read.count.clone.2,
                "el4.rank" = el4.rank, "el4.count" = el4.count, "el4.reads" = read.count.clone.3)
    contamination.qc[i,] <- qc.row
  
    ##
    ## Recalculation - re-calculate clone fraction.
    ##
  
    new_total_count_v <- sum(new.clone[[rawCountCol_v]])
    new.clone[[fracCol_v]] <- new.clone[[rawCountCol_v]] / new_total_count_v
    
    ## Count clean clone files as well
    clean.counter <- 0
    if (nrow(new.clone) == nrow(curr.clone)){
        clean.counter <- clean.counter + 1
        cat("No contaminant sequences found in sample:", i, "\n")
    }
    
    ## Remove extension from file name and add "_decontam"
    out.name <- paste(gsub(".txt", '', clone.files[i]), "_decontam.txt", sep = '')
  
    ## Write output to new file
    if (which.outputs %in% c("clones", "both")){
        ## Write output
        write.table(new.clone, paste(out.dir, out.name, sep = ''),
                    sep = '\t', quote = F, row.names = F)
        ## Update progress
        print(c("Writing output to ", paste(out.dir, out.name, sep = '')))
    } # fi

    ## Update progress
    if (i %% 10 == 0){print(c("Currently on clone ", i))}
    
} # for

colnames(contamination.qc) <- c("Sample", "Orig.Unique.Clones", "Orig.Total.Clones", "Remaining.Clones", "Contam.clones",
                                "p14.rank", "p14.count", "p14.reads", "ot1.rank", "ot1.count", "ot1.reads", "el4.rank", "el4.count", "el4.reads")

cat(clean.counter, "clone files were uncontaminated!\n")


if (which.outputs %in% c("qc", "both")){
    ## Write QC output
    write.table(contamination.qc, file = paste0(qc.dir, batch_v, "_contaminationQC.txt"),
                sep = '\t', quote = F, row.names = F)
    ## Update
    print(c("Writing QC output to ", paste0(qc.dir, batch_v, "_contaminationQC.txt")))
} # fi
