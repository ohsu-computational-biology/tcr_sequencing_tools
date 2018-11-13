#!/bin/bash

# Run all of the QC scripts at one time

SCRIPTS=$tool/QC/
REF=$tool/reference/
DATA=$data/
OUT=$data/QC/

echoerr() { printf "%s\n" "$*" >&2; }

# ALIGN
echoerr align
Rscript $SCRIPTS/mixcr.alignment.QC.R $DATA/mixcr/reports/align/ $OUT/

# ASSEMBLE
echoerr assemble
Rscript $SCRIPTS/mixcr.assemble.QC.R $DATA/mixcr/reports/assemble/ $OUT/

