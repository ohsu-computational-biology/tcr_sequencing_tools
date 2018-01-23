#!/bin/sh

### Organize TCRseq fastq files into multiple sub-directories for easier import into Galaxy

DIR=$1 # directory of fastq files in dataset_library_import directory. Must be batch name, i.e. first portion of all the file names. Must have trailing slash.
BATCH=`basename $DIR`
BASE=${DIR%$BATCH/}

echo "DIR: " $DIR
echo "BATCH: " $BATCH
echo "BASE: " $BASE

### Determine number of subdirectories to make
len=`ls $DIR | wc -l`
len=`expr $len / 20`


for i in $(seq 0 $len); do

    ## Get second part of directory
    sec=`expr $i + 1` 

    ## Skip if odd; move files if even
    if [ `expr $i % 2` == 1 ]; then
	continue
    else
	## Make directory
	CURRDIR=$DIR/$i\_$sec
	mkdir -p $CURRDIR
	echo "Currently on: " $i
	echo "Directory: " $CURRDIR
	## Move files there
	for j in {0..9}; do
	    ## Move first half. Different if first one, because no tens place
	    if [ $i == 0 ]; then
	        mv $DIR/$BATCH\_S$j\_R* $CURRDIR
	    else
		mv $DIR/$BATCH\_S$i$j\_R* $CURRDIR
	    fi
	    ## Move second half.
	    mv $DIR/$BATCH\_S$sec$j\_R* $CURRDIR
	done
	## Move to their own directory within your dataset_library_import, so galaxy can see them.
	NEWDIR=$BASE/$BATCH\_$i\_$sec
	echo "Moving " $CURRDIR " to " $NEWDIR
	mv $CURRDIR $NEWDIR
    fi
done	


