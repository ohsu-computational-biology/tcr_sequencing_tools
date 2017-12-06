#!/bin/sh

### Create symlinks for all of the data in a particular directory to another directory.
### Intended use is to link all of the fastq or fastq.gz files for a particular TCRseq batch
### into a galaxy import directory.



### VARIABLES
TCR=/home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/
GALAXY=/home/exacloud/lustre1/Galaxy/dataset_library_import/

### ARGUMENTS
IN=$1		# Sub-directory within $TCR that contains fastq files. Usually DNAXXXXLC/fastqs_from_core/fastqs or similar
OUT=$2          # Sub-directory within $GALAXY to link files to. Usually user@ohsu.edu/DNAXXXXLC where DNAXXXXLC is same as in $IN

### BODY

cd $GALAXY/$OUT/

for file in `ls $TCR/$IN`; do
    ln -s $TCR/$IN/$file $file
done

