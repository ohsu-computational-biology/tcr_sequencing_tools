#!/bin/sh

DIR=$1
rename_log=$2

### Change spike count file from csv to tsv
for file in `ls $DIR`; do
    cat $DIR/$file | tr "," "\t" > temp;
    mv -f temp $DIR/$file;
done

### Remove extra column of row numbers
for file in `ls $DIR`; do
    cat $DIR/$file | awk -F '\t' -v OFS='\t' '{if (NR == 1) print $0; else print $2, $3, $4, $5, $6}' > temp;
    mv -f temp $DIR/$file;
done

### Remove extra sample number in file names
if [[ $rename_log == "TRUE" ]]; then
   for file in `ls $DIR`; do
       new=`echo $file | sed -E 's/(^.*_)(S[0-9]+_)(S[0-9]+)(\..*$)/\1\3\4/'`;
       mv -f $DIR/$file $DIR/$new;
   done
fi   
