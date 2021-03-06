###################################################
### Calculate Negative-Binomial Scaling Factors ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################################################

### Summary:
  ### Use data from spike-only samples to calculate more robust scaling factors that can be applied to the raw clonotype counts.
  ### This particular instance has 3 sets of spike-only samples. For each set of samples, read in the 25-bp spike counts,
  ### then calculate the dispersion (theta) and estimate (mu) of the negative binomial model that is fit to the counts. Take each
  ### set of estimates (260 - one for each VJ combination) and combine them into 1 set of estimates. Then take the mean of that set.
  ### Finally, divide each combined estimate by the mean in order to get the final scaling factor. Output is a data.frame containing
  ### the 260 final scaling factors written to a file.

### Standard Formula: normCount(i) = rawCount(i) * (mean_mu /  mu(i))
  ### where mean_mu is the mean of the mean of all of the dispersion estimates
    ### We take 3 sets of spike-only samples and estimate dispersions for each (3 sets of 260 dispersions)
    ### Then we take the mean of those (1 set of 260 dispersions), then we take the mean of that (1 value).
    ### mu(i) then is one of the 260 mu's from the 1 set of 260 dispersions.

### Note:
  ### The scaling factor calculated here is mu(i) / mean_mu for each mu. As shown above, the actual formula is mean_mu / mu(i) for
  ### each mu. What this means is that in order to adjust a raw count to a normalized count using the values output by this script,
### one must use the formula: normCount(i) = rawCount(i) / scalingFactor(i); where scalingFactor(i) = mu(i) / mean_mu

### Final note:
### Changed the scaling factor calculation to mean_mu / mu(i) for each mu, so these scaling factors act the same as the
### naive normalization method!

####################
### Dependencies ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################
library(data.table)
library(MASS)
library(plyr)

#################
### Arguments ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################

arguments <- commandArgs(trailingOnly = T)
dataDir_v <- arguments[1] # base directory that contains folders w/ 25-bp spike counts
outDir_v <- arguments[2] # directory to write out scaling factor file
check_v <- arguments[3] # T or F. If T, plot adjusted counts of one of the inputs.

### For testing
# dataDir_v <- "/Users/hortowe/Desktop/OHSU/tcr_spike/data/nb_norm/"
# outDir_v <- "/Users/hortowe/Desktop/OHSU/tcr_spike/data/nb_norm/output/"
# check_v <- F

########################
### Define Functions ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########################

### Read in spike data
readSpikes <- function(file_v,directory_v, sep_v = "\t") {
    
  ## This function is used within readDir via apply()  
  ## file_v - specific file name within a data directory. 
  ## directory_v - directory containing spike count files
  ## sep_v - column separator of file_v. Usually \t, but can be , sometimes
    
  ## Get spike data
  spike_dt <- fread(file.path(directory_v, file_v), sep = sep_v)
  ## Subset to contain only V, J, and spike counts
  spike_dt <- spike_dt[,c(3:5), with = F]
  ## Subset name to only include sample number
  sample_v <- (gsub(".*_S|*.assembled.*", '', file_v))
  ## Add column names
  colnames(spike_dt) <- c("V","J",sample_v)
  ## Fix column data (trailing - after V genes)
  spike_dt[,V := gsub('-', '', V)]
  ## Return spike
  return(spike_dt)
} # readSpikes

### Read in directory of data
readDir <- function(sub_v) {

  ## sub_v - sub-directory of the global variable dataDir_v
    
  ## Get directory
  directory_v <- paste0(dataDir_v, sub_v)
  ## Get files in directory
  allFiles_v <- list.files(directory_v, pattern = "txt$")
  ## Read in files using readSpikes
  spikes_dt <- lapply(allFiles_v, readSpikes, directory_v = directory_v)
  ## Merge together
  spikes_dt <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("V", "J"), all = T), spikes_dt)
  ## Convert V and J columns to rownames
  spikes_dt <- colToRowNames(spikes_dt)
  ## Sort
  setcolorder(spikes_dt, as.character(sort(as.numeric(colnames(spikes_dt)))))
  ## Return
  return(spikes_dt)
} # readDir

### Change V and J columns to rownames
colToRowNames<-function(data_dt, cols_v = c("V", "J")) {

  ## data_dt - a data.table object of spike counts
  ## cols_v - column names to be used to extract values for new  rownames  

  ## Paste selected columns together and assign as row names
  rownames(data_dt) <- paste0(data_dt[[cols_v[1]]], " | ", data_dt[[cols_v[2]]])
  ## Remove those columns
  data_dt[,(cols_v) := NULL]
  ## Convert NAs to 0
  data_dt[is.na(data_dt)] <- 0
  ## Return
  return(data_dt)
}

### Calcualte dispersion
calcDispersion <- function(spike_v, cols_v){

  ## spike_v - vector of 260 spike counts
  ## cols_v - vector of columns to extract from NB model
    
  ## Create Negative Binomial Model
  m.nb <- glm.nb (spike_v ~ 1)
  ## Extract theta
  theta_v <- summary(m.nb)$theta
  ## Extract mu
  mu.nb <- coef(summary(m.nb, dispersion = 1))
  ## Calculate results
  results_dt <- data.frame(mu.nb, theta_v)
  ## Paste batch prefix to names, except entirely new name for theta
  colnames(results_dt) <- c(paste(cols_v[1], colnames(results_dt)[1:4], sep = '.'), cols_v[2])
  ## Return
  return(results_dt)
} # calcDispersion

################
### Get Data ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################

### Read in counts
### TODO - this is hardcoded...
newCount_dt <- readDir(sub_v = "counts/LIB170111LC/")
oldCount_dt <- readDir(sub_v = "counts/DNA160609LC/")

### Subset data
oldCount_df <- as.data.frame(oldCount_dt[,1:20]); rownames(oldCount_df) <- rownames(oldCount_dt)
newDil_df <- as.data.frame(newCount_dt[,1:10]); rownames(newDil_df) <- rownames(newCount_dt)
dil5_df <- as.data.frame(newCount_dt[,11:20]); rownames(dil5_df) <- rownames(newCount_dt)

######################################################
### Calculate Dispersion (theta) and Estimate (mu) ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################################################

### Calculate dispersion for each sample, combine to data.frame, remove extraneous column
oldDisp_lsdf <- apply(oldCount_df, 1, function(x) calcDispersion(unlist(x), cols_v = c("mu.nb_old", "theta_old")))
oldDisp_df <- ldply(oldDisp_lsdf, rbind); oldDispNB_df <- oldDisp_df[,-1]

newDilDisp_lsdf <- apply(newDil_df, 1, function(x) calcDispersion(unlist(x), cols_v = c("mu.nb_newdil", "theta_newdil")))
newDilDisp_df <- ldply(newDilDisp_lsdf, rbind); newDilDispNB_df <- newDilDisp_df[,-1]

dil5Disp_lsdf <- apply(dil5_df, 1, function(x) calcDispersion(unlist(x), cols_v = c("mu.nb_dil5", "theta_dil5")))
dil5Disp_df <- ldply(dil5Disp_lsdf, rbind); dil5DispNB_df <- dil5Disp_df[,-1]

### Combine all of the data
### TODO - is this necessary?
allCount_dt <- cbind(newCount_dt, newDilDispNB_df, dil5DispNB_df, oldCount_dt, oldDispNB_df)

### Extract mu (estimate) and theta (dispersion)
allEstimates_dt <- allCount_dt[,c("mu.nb_newdil.Estimate", "mu.nb_dil5.Estimate", "mu.nb_old.Estimate",
                                  "theta_old", "theta_newdil", "theta_dil5")]

### Calculate estimate differences and combine them
newDilMinusOldEstimate_v <- (sum(allEstimates_dt$mu.nb_newdil.Estimate) - sum(allEstimates_dt$mu.nb_old.Estimate)) / 260
dil5MinusOldEstimate_v <- (sum(allEstimates_dt$mu.nb_dil5.Estimate) - sum(allEstimates_dt$mu.nb_old.Estimate)) / 260
combinedEstimate_v <- ((2*allEstimates_dt$mu.nb_old.Estimate) +                               # Weight old estimate b/c 20?
                       allEstimates_dt$mu.nb_newdil.Estimate - newDilMinusOldEstimate_v +   # newDil estimate - diff b/w newDil and old
                       allEstimates_dt$mu.nb_dil5.Estimate - dil5MinusOldEstimate_v) / 4    # dil5 estimate - diff b/w dil5 and old
                                                                                            # divide by 4 to take the mean?
allEstimates_dt$"mu.nb_combined.Estimate" <- combinedEstimate_v

#################################
### Scaling Factor and Output ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################################

### Turn estimates into data.frame w/ appropriate VJ row names
combinedEstimate_df <- data.frame(combinedEstimate_v); rownames(combinedEstimate_df) <- rownames(oldCount_df)

### Apply exponential function to counts (new.count = e^old.count)
combinedEstimate_df <- exp(combinedEstimate_df)

### Take arithmetic mean
meanEstimate_v <- colMeans(combinedEstimate_df)

### "Normalize" estimates by dividing by mean
                                        #normEstimate_df <- combinedEstimate_df / meanEstimate_v
normEstimate_df <- meanEstimate_v / combinedEstimate_df

### Add extra set for V121
rownames(normEstimate_df) <- gsub("V1212", "V121", rownames(normEstimate_df))
temp122 <- normEstimate_df[grep("V121", rownames(normEstimate_df)), ,drop = F ]
rownames(temp122) <- gsub("V121", "V122", rownames(temp122))
normEstimate_df <- rbind(normEstimate_df, temp122)

### Add V and J columns (for compatibility with other normalization method)
outputV_v <- gsub(" \\| J.*$", "", rownames(normEstimate_df))
outputJ_v <- gsub("^.* \\| ", "", rownames(normEstimate_df))
output_df <- cbind(outputV_v, outputJ_v, normEstimate_df)
colnames(output_df) <- c("V", "J", "scaling.factor")

### Write out scaling factors
write.table(output_df, file = file.path(outDir_v, "nb.scaling.factors.txt"), 
            sep = '\t', quote = F, row.names = F)

##############
### Check! ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############

if (check_v == T){
  ### Convert to data.table and change rownames to column
  setDT(oldCount_df, keep.rownames = T); colnames(oldCount_df)[1] <- "VJ"
  
  ### Add normalized estimates (scaling factors)
  oldCount_df$mu <- normEstimate_df
  
  ### Melt for plotting
  oldCountMelt_df <- melt(oldCount_df, id.vars = c("VJ", "mu"))
  
  ### Adjust counts
  oldCountMelt_df$"adj.Value" <- oldCountMelt_df$value / oldCountMelt_df$mu
  
  ### Rename columns and remove estimate
  colnames(oldCountMelt_df)[c(3,4)] <- c("sample", "unnorm")
  oldCountMelt_df <- oldCountMelt_df[,c("VJ", "sample", "unnorm", "adj.Value")]
  
  ### Convert to double
  oldCountMelt_df$unnorm <- as.double(oldCountMelt_df$unnorm)
} # fi
