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

### Sort files
clone.files <- list.files(clone.dir)
clone.files <- clone.files[order(as.numeric(gsub(".*_S|_align.*", '', clone.files)))]

### Empty matrix for output summary
contam_reads_mat <- matrix(nrow = length(clone.files), ncol = 2)

### Batch name for output files
batch_v <- unlist(strsplit(clone.files[1], split = "_"))[1]

### Empty QC Matrix
contamination.qc <- matrix(nrow = length(clone.files), ncol = 11)

###################
### Calculation ###
###################

# Iterate through each clone file, removing monoclonal contaminants
for (i in 1:length(clone.files)){
    ## Counter for how many reads are devoted to contaminated sequences
    curr_contam_count_v <- 0
    ## Sample name
    currName_v <- unlist(strsplit(clone.files[i], split = "_"))[2]
    
    ## Get a clone file
    curr.clone <- fread(paste(clone.dir, clone.files[i], sep = ''))

    ## Add index for ranking
    curr.clone$index <- seq(1, length(curr.clone$`Clone ID`))

    ## Get first QC info
    orig.unique.count <- length(curr.clone$`Clone ID`)
    orig.total.count <- sum(curr.clone$`Clone count`)
    

    ## Update column names and contents for V and J segments
    curr.clone$`V segments` <- gsub("TRB|\\*00", "", curr.clone$`Best V hit`)
    curr.clone$`V segments` <- gsub("-", "", curr.clone$`V segments`)
    curr.clone$`J segments` <- gsub("TRB|\\*00", "", curr.clone$`Best J hit`)

    
    ## Extract clones (take just the first entry if more than one clone matches the criteria)
    ## p14
    offending.clone.1 <- curr.clone[(curr.clone$`V segments` == "V133" # p14
        & curr.clone$`J segments` == "J2-4" &
          curr.clone$`AA. Seq. CDR3` == "CASSDAGGRNTLYF"),]
    
    offending.clone.1 <- offending.clone.1[1,] # subset for first only
    
    ## Update count of contaminated sequences
    if (is.na(offending.clone.1$Reads[1])) {
        curr_contam_count_v <- curr_contam_count_v + 0
    } else {
        curr_contam_count_v <- curr_contam_count_v + length(unlist(strsplit(offending.clone.1$Reads, split = ',')))
    }


    ## ot1
    offending.clone.2 <- curr.clone[(curr.clone$`V segments` == "V121" # ot1
        & curr.clone$`J segments` == "J2-7" &
          curr.clone$`AA. Seq. CDR3` == "CASSRANYEQYF"),]
    offending.clone.2 <- offending.clone.2[1,] # subset for first only

    ## Update count of contaminated sequences
    if (is.na(offending.clone.2$Reads[1])) {
        curr_contam_count_v <- curr_contam_count_v + 0
    } else {
        curr_contam_count_v <- curr_contam_count_v + length(unlist(strsplit(offending.clone.2$Reads, split = ',')))
    }


    ## el4
    offending.clone.3 <- curr.clone[(curr.clone$`V segments` == "V15" # el4
        & curr.clone$`J segments` == "J2-3" &
          curr.clone$`AA. Seq. CDR3` == "CASSTGTETLYF"),]
    offending.clone.3 <- offending.clone.3[1,] # subset for first only
  
    if (is.na(offending.clone.3$Reads[1])) {
        curr_contam_count_v <- curr_contam_count_v + 0
    } else {
        curr_contam_count_v <- curr_contam_count_v + length(unlist(strsplit(offending.clone.3$Reads, split = ',')))
    }

    ## Add count of contaminated reads to output matrix
    contam_reads_mat[i,] <- c(currName_v, curr_contam_count_v)

    
    ## Combine contaminants together into a data frame with the original clone data frame
    combo <- rbind(curr.clone, offending.clone.1, offending.clone.2, offending.clone.3)
    ## Remove duplicated entries (i.e. offending clones) from combined dataframe
    new.clone <- combo[!duplicated(combo, fromLast = T) & seq(nrow(combo)) <= nrow(curr.clone),]
    new.clone <- new.clone[!(is.na(new.clone$`Clone ID`)),]

    ## QC Data
    p14.rank <- offending.clone.1$index[1]
    p14.count <- sum(offending.clone.1$`Clone count`)
    ot1.rank <- offending.clone.2$index[1]
    ot1.count <- sum(offending.clone.2$`Clone count`)
    el4.rank <- offending.clone.3$index[1]
    el4.count <- sum(offending.clone.3$`Clone count`)
    remaining.count <- orig.total.count - p14.count - ot1.count - el4.count
    qc.row <- c("Sample" = currName_v, "Orig.Unique.Clones" = orig.unique.count, "Orig.total.Clones" = orig.total.count,
                "Remaining.Clones" = remaining.count, "Contam.clones" = curr_contam_count_v, "p14.rank" = p14.rank,
                "p14.count" = p14.count, "ot1.rank" = ot1.rank, "ot1.count" = ot1.count, "el4.rank" = el4.rank, "el4.count" = el4.count)
    contamination.qc[i,] <- qc.row
  
    ##
    ## Recalculation - re-calculate clone fraction.
    ##
  
    new_total_count_v <- sum(new.clone$`Clone count`)
    new.clone$`Clone fraction` <- new.clone$`Clone count` / new_total_count_v
    
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

cat(clean.counter, "clone files were uncontaminated!\n")

write.table(contam_reads_mat, paste0(qc.dir, batch_v, "_contam_clone_reads.txt"),
            sep = '\t', quote = F, row.names = F)


if (which.outputs %in% c("qc", "both")){
    ## Write QC output
    write.table(contamination.qc, file = paste0(qc.dir, batch_v, "_contaminationQC.txt"),
                sep = '\t', quote = F, row.names = F)
    ## Update
    print(c("Writing QC output to ", paste0(qc.dir, batch_v, "_contaminationQC.txt")))
} # fi
