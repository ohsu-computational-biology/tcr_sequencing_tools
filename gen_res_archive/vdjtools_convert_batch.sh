#!/bin/bash

WORKING_DIR="/mnt/lustre1/CompBio/data/tcrseq/dhaarini/DNA150826LC/150924_NS500681_0014_AH3NTVAFXX/postprocessing/vdjtools_converted"

echo "Running VDJTools Convert() in batch form in" $WORKING_DIR

for i in {1..1}; do
	CUR_FILE=`ls $WORKING_DIR | grep DNA150826LC_${i}_`
	echo "Current file: " $CUR_FILE
	java -Xmx4G -jar /mnt/lustre1/CompBio/genomic_resources/tcrseq/vdjtools/vdjtools-1.0.2.jar Convert \
        -S mixcr \
        $WORKING_DIR/$CUR_FILE \
        $WORKING_DIR
done
