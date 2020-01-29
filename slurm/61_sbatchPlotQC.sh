#!/bin/sh

#SBATCH --partition          exacloud                # partition (queue)
#SBATCH --nodes              1                       # number of nodes
#SBATCH --ntasks             1                       # number of "tasks" to be allocated for the job
#SBATCH --ntasks-per-core    1                       # Max number of "tasks" per core.
#SBATCH --cpus-per-task      1                       # Set if you know a task requires multiple processors
##SBATCH --mem-per-cpu        8000                    # Memory required per allocated CPU (mutually exclusive with mem)
#SBATCH --mem                16000                  # memory pool for each node
#SBATCH --time               0-24:00                 # time (D-HH:MM)
#SBATCH --output             newQC_%j.out        # Standard output
#SBATCH --error              newQC_%j.err        # Standard error

### Directories
TOOL=$tool/51_plotQC
DATA=$data
PROJ=`basename $DATA`
OUT=$DATA/QC/plots

mkdir -p $OUT

### Toggle
PEAR=true
CLONE=true
CONTAM=true
SPIKE=true
FASTQC=true

### PEAR
if $PEAR; then
	Rscript $TOOL/pearQC.R \
		--inputDir $DATA/QC/pear/ \
		--outDir $OUT
fi

### READS AND CLONES
if $CLONE; then
	Rscript $TOOL/readAndCloneQC.R \
		--contamFile $DATA/QC/std/$PROJ\_contaminationQC.txt \
		--nineFile $DATA/QC/std/count.spikes.9bp.QC.summary.txt \
		--outDir $OUT \
		--plot 'dot'
fi

### READS AND CONTAM
if $CONTAM; then
	Rscript $TOOL/readAndContamQC.R \
		--contamFile $DATA/QC/std/$PROJ\_contaminationQC.txt \
		--nineFile $DATA/QC/std/count.spikes.9bp.QC.summary.txt \
		--outDir $OUT
fi

### SPIKE
if $SPIKE; then
	Rscript $TOOL/spikeQC.R \
		--nine $DATA/spike_counts/9bp/qc \
		--twentyFive $DATA/spike_counts/25bp/qc \
		--outDir $OUT
fi

### FASTQC
if $FASTQC; then
	Rscript $TOOL/fastQC_sequenceQuality.R \
	--inputFile $DATA/fastqs_from_core/FastQC/multiqc_data/multiqc_fastqc.txt \
	--outDir $OUT
fi

