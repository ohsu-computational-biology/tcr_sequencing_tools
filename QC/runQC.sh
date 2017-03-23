#!/bin/bash

# Run all of the QC scripts at one time

SCRIPTS=$tool/QC/
REF=$tool/reference/
DATA=$data/
OUT=$data/QC/

# REMOVE SPIKES
echo "Remove spikes"
Rscript $SCRIPTS/remove.spikes.QC.R $DATA/peared_fastqs/assembled/ $DATA/mixcr/despiked_fastqs/ $OUT/

# COUNT SPIKES
echo "Count spikes"
Rscript $SCRIPTS/count.spikes.QC.R $DATA/spike_counts/9bp/qc/ $OUT/
mv $OUT/count.spikes.QC.summary.txt $OUT/count.spikes.9bp.QC.summary.txt

Rscript $SCRIPTS/count.spikes.QC.R $DATA/spike_counts/25bp/qc/ $OUT/
mv $OUT/count.spikes.QC.summary.txt $OUT/count.spikes.25bp.QC.summary.txt

# ALIGN
echo "align"
Rscript $SCRIPTS/mixcr.alignment.QC.R $DATA/mixcr/reports/align/ $OUT/

# ASSEMBLE
echo "assemble"
Rscript $SCRIPTS/mixcr.assemble.QC.R $DATA/mixcr/reports/assemble/ $OUT/

# NON-STANDARD ALIGNMENTS
echo "non-standard alignments"
Rscript $SCRIPTS/countNonStandardHits.R $DATA/normalization/decontam $REF/text_barcodesvj.txt $OUT/

# NORMALIZATION
echo "normalization"
Rscript $SCRIPTS/normalization.QC.R $DATA/normalization/decontam/ $DATA/normalization/normalized_clones/ $DATA/normalization/QC/

Rscript $SCRIPTS/aggregate.normalization.QC.R $DATA/normalization/QC/ $OUT/
