####################################################
### Remove monoclonal contamination from samples ###
####################################################

library(data.table)


##############
### Set Up ###
##############

arguments <- commandArgs(trailingOnly = T)

clone.dir <- arguments[1]        # #data/mixcr/export_clones
out.dir <- arguments[2]

### Sort files
clone.files <- list.files(clone.dir)
clone.files <- clone.files[order(as.numeric(gsub(".*_S|_align.*", '', clone.files)))]


###################
### Calculation ###
###################

# Iterate through each clone file, removing monoclonal contaminants
for (i in 1:length(clone.files)){
    ## Get a clone file
    curr.clone <- fread(paste(clone.dir, clone.files[i], sep = ''))

    ## Update column names and contents for V and J segments
    curr.clone$`V segments` <- gsub("TRB|\\*00", "", curr.clone$`Best V hit`)
    curr.clone$`V segments` <- gsub("-", "", curr.clone$`V segments`)
    curr.clone$`J segments` <- gsub("TRB|\\*00", "", curr.clone$`Best J hit`)
    
    ## Remove extension from file name and add "_decontam"
    out.name <- paste(gsub(".txt", '', clone.files[i]), "_new_cols.txt", sep = '')
  
    ## Write output to new file
    write.table(curr.clone, paste(out.dir, out.name, sep = ''),
                sep = '\t', quote = F, row.names = F)
    ## Update progress
    print(c("Writing output to ", paste(out.dir, out.name, sep = '')))

} # for
