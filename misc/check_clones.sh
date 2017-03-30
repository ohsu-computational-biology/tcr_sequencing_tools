#!/bin/sh

DIR=$1

for file in $DIR/*; do
    start_cols=`eval head -1 $file | sed 's/\t/\n/g' | wc -l`
    end_cols=`eval tail -1 $file | sed 's/\t/\n/g' | wc -l`
    col_diff=`expr $start_cols - $end_cols`
    if [ $col_diff != 0 ]; then
	echo $file "is truncated. First column has: " $start_cols " while last column has: " $end_cols
    fi
done
