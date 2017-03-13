#!/bin/sh

# This script takes replicates that are offset by a specific amount and combines them together

firstSet=$1
offset=$2

# Set Up
cd $data
mkdir -p ./peared_fastqs/combine_rep

temp_file=`eval ls ./peared_fastqs/assembled/ | head -1`
batch=${temp_file%\_*}

# Run Combination
cd $data/peared_fastqs/assembled

for i in $(seq 1 $firstSet); do
    
    # Get replicate number
    rep=$(($i + $offset))

    # Output name
    out_file=$batch\_S$i\.assembled.combined.S$rep\.fastq
    
    # Copy first replicate
    eval `cp ./$batch\_S$i\.assembled.fastq ../combine_rep/$out_file`

    # Add second replicate
    eval `cat ./$batch\_S$rep\.assembled.fastq >> ../combine_rep/$out_file`

    # Print numbers
    echo $i $rep $out_file
    printf "First: %s\tSecond: %s\tOut: %s\n" "$i" "$rep" "$out_file" >> $data/QC/rep.combine.log

done
