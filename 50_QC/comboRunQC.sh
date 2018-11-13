#!/bin/bash

# Run all of the QC scripts at one time

SCRIPTS=$tool/QC/
DATA=$data/
OUT=$data/combo_QC/

# REMOVE SPIKES
echo "Remove spikes"
Rscript $SCRIPTS/remove.spikes.QC.R $DATA/peared_fastqs/combine_rep/ $DATA/combo_mixcr/despiked_fastqs/ $OUT/

# COUNT SPIKES
echo "Count spikes"
#Rscript $SCRIPTS/count.spikes.QC.R $DATA/spike_counts/combo/9bp/qc/ $OUT/
#mv $OUT/count.spikes.QC.summary.txt $OUT/count.spikes.9bp.QC.summary.txt

Rscript $SCRIPTS/count.spikes.QC.R $DATA/spike_counts/combo/25bp/qc/ $OUT/
mv $OUT/count.spikes.QC.summary.txt $OUT/count.spikes.25bp.QC.summary.txt

# ALIGN
echo "align"
#Rscript $SCRIPTS/mixcr.alignment.QC.R $DATA/combo_mixcr/reports/align/ $OUT/

# ASSEMBLE
echo "assemble"
#Rscript $SCRIPTS/mixcr.assemble.QC.R $DATA/combo_mixcr/reports/assemble/ $OUT/

# NORMALIZATION
echo "normalization"
#Rscript $SCRIPTS/normalization.QC.R $DATA/combo_normalization/decontam/ $DATA/combo_normalization/normalized_clones/ $DATA/combo_normalization/QC/

#Rscript $SCRIPTS/aggregate.normalization.QC.R $DATA/combo_normalization/QC/ $OUT/
