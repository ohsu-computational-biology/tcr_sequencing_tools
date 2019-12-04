#!/bin/bash

# Run all of the QC scripts at one time

SCRIPTS=$tool/50_QC/
REF=$tool/00_reference/
DATA=$data/
OUT=$data/QC/

echoerr() { printf "%s\n" "$*" >&2; }

# REMOVE SPIKES
echoerr Remove spikes
/usr/bin/Rscript $SCRIPTS/remove.spikes.QC.R $DATA/peared_fastqs/assembled/ $DATA/mixcr/despiked_fastqs/ $OUT/

# COUNT SPIKES
echoerr Count spikes 9
Rscript $SCRIPTS/count.spikes.QC.R $DATA/spike_counts/9bp/qc/ $OUT/
mv $OUT/count.spikes.QC.summary.txt $OUT/count.spikes.9bp.QC.summary.txt

echoerr Count spikes 25
Rscript $SCRIPTS/count.spikes.QC.R $DATA/spike_counts/25bp/qc/ $OUT/
mv $OUT/count.spikes.QC.summary.txt $OUT/count.spikes.25bp.QC.summary.txt

# MIXCR
echoerr analyze
Rscript $SCRIPTS/mixcr.analyze.QC.R --inputDir $DATA/mixcr/reports/ --outDir $OUT/

# NON-STANDARD ALIGNMENTS
echoerr non-standard alignments
Rscript $SCRIPTS/countNonStandardHits.R $DATA/normalization/decontam/ $REF/text_barcodesvj.txt $OUT/

# NORMALIZATION
echoerr normalization
Rscript $SCRIPTS/normalization.QC.R -r $DATA/normalization/decontam/ -n $DATA/normalization/normalized_clones/ -o $DATA/normalization/QC/

