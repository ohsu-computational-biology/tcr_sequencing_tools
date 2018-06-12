###########################################################
### Add new columns V and J segments without extra bits ###
###########################################################

library(data.table)


##############
### Set Up ###
##############

arguments <- commandArgs(trailingOnly = T)

clone.dir <- arguments[1]        # #data/mixcr/export_clones
out.dir <- arguments[2]
rename_v <- arguments[3]

### Sort files
clone.files <- list.files(clone.dir)
clone.files <- clone.files[order(as.numeric(gsub(".*_S|_align.*|_clones.*", '', clone.files)))]


###################
### Calculation ###
###################

# Iterate through each clone file, removing monoclonal contaminants
for (i in 1:length(clone.files)){
    ## Get a clone file
    curr.clone <- fread(file.path(clone.dir, clone.files[i]))

    ## Get columns
    vCol_v <- grep('[Bb]est[ ]*V[ ]*[Hh]it', colnames(curr.clone), value = T)
    jCol_v <- grep('[Bb]est[ ]*J[ ]*[Hh]it', colnames(curr.clone), value = T)

    ## Update column names and contents for V and J segments
    curr.clone$`V segments` <- gsub("TRB|\\*00", "", curr.clone[[vCol_v]])
    curr.clone$`V segments` <- gsub("-", "", curr.clone$`V segments`)
    curr.clone$`J segments` <- gsub("TRB|\\*00", "", curr.clone[[jCol_v]])
    
    ## Remove extension from file name and add "_decontam"
    if (rename_v) {
        out.name <- paste(gsub(".txt", '', clone.files[i]), "_new_cols.txt", sep = '')
    } else {
	out.name <- clone.files[i]
    }
 
    ## Write output to new file
    write.table(curr.clone, file.path(out.dir, out.name),
                sep = '\t', quote = F, row.names = F)
    ## Update progress
    print(c("Writing output to ", paste(out.dir, out.name, sep = '')))

} # for
