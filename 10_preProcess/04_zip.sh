#!/bin/bash

#PROJ=$data
#INDIR=$PROJ/fastqs_from_core/fastqs

#echo $INDIR

SAMPLE=$1
#FILE=$INDIR/$SAMPLE

gzip $SAMPLE

