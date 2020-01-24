#!/bin/sh

### Directories
TOOL=/Users/hortowe/my_tool_repos/tcr_sequencing_tools/newQC
DATA=/Users/hortowe/OHSU/tcr_spike/data/LIB191025LC
PROJ=LIB191025LC

### Toggle
PEAR=false
CLONE=false
CONTAM=false
SPIKE=true
FASTQC=false

### PEAR
if $PEAR; then
	Rscript $TOOL/pearQC.R \
		--inputDir $DATA/qc_inputs/pear/ \
		--outDir $DATA/qc_outputs/pear/
fi

### READS AND CLONES
if $CLONE; then
	Rscript $TOOL/readAndCloneQC.R \
		--contamFile $DATA/qc_inputs/reads_and_contam/$PROJ\_contaminationQC.txt \
		--nineFile $DATA/qc_inputs/reads_and_contam/count.spikes.9bp.QC.summary.txt \
		--outDir $DATA/qc_outputs/reads_and_contam \
		--plot 'dot'
fi

### READS AND CONTAM
if $CONTAM; then
	Rscript $TOOL/readAndContamQC.R \
		--contamFile $DATA/qc_inputs/reads_and_contam/$PROJ\_contaminationQC.txt \
		--nineFile $DATA/qc_inputs/reads_and_contam/count.spikes.9bp.QC.summary.txt \
		--outDir $DATA/qc_outputs/reads_and_contam
fi

### SPIKE
if $SPIKE; then
	Rscript $TOOL/spikeQC.R \
		--nine $DATA/qc_inputs/spike9 \
		--twentyFive $DATA/qc_inputs/spike25 \
		--outDir $DATA/qc_outputs/spike
fi

### FASTQC
if $FASTQC; then
	Rscript $TOOL/fastQC_sequenceQuality.R \
	--inputFile $DATA/FastQC/multiqc_data/multiqc_fastqc.txt \
	--outDir $DATA/qc_outputs/fastqc
fi

