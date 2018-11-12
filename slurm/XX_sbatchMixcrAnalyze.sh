#!/bin/bash

# This is new for MiXCR 3.0 - the `analyze` module will run all 3 of `align`, `assemble`, and `exportClones` at once.

#SBATCH --partition          exacloud                   # partition (queue)
#SBATCH --nodes              1                          # number of nodes
#SBATCH --ntasks             1                          # number of "tasks" to be allocated for the job
#SBATCH --ntasks-per-core    1                          # Max number of "tasks" per core.
#SBATCH --cpus-per-task      1                          # Set if you know a task requires multiple processors
#SBATCH --mem                32000                      # memory pool for each node
#SBATCH --time               0-24:00                    # time (D-HH:MM)
#SBATCH --output             mixcr_analyze_%A_%a.out    # Standard output
#SBATCH --error              mixcr_analyze_%A_%a.err    # Standard error
#SBATCH --array              1-1                        # sets number of jobs in array

### Set I/O variables

IN=$data/mixcr/despiked_fastqs             # Directory containing all input files. Should be one job per file
OUT=$data/mixcr/           # Directory where output files should be written
REPORTDIR=$data/mixcr/reports/
INDEX=$data/mixcr/indexes
MYBIN=$MIXCR/mixcr.jar          # Path to shell script or command-line executable that will be used
PRESET=$tool/mixcr/analyze_preset.txt
SPECIES="mmu"
TEST=False

### Record slurm info

date
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
printf "\n\n"

#################
### VARIABLES ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################

### todo file
TODO=$data/tools/todo/analyze.txt
CURRFILE=`awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}' $TODO`
BASE="${CURRFILE%%S[0-9]*}"
SNUM="${CURRFILE%%.assembled.removed.fastq}"; SNUM="${SNUM##*S}"

### Issign I/O files/paths
FASTQ=$IN/$CURRFILE
REPORT=S$SNUM\_report.txt
EXPORT=$BASE\S$SNUM

###############
### ANALYZE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###############

### Run MiXCR
mixcr analyze amplicon \
	-s mmu \
	--starting-material dna \
	--5-end v-primers \
	--3-end j-primers \
	--adapters no-adapters \
	--report $REPORTDIR/$REPORT \
	--receptor-type trb \
	--region-of-interest CDR3 \
	--only-productive  \
	--align "-OsaveOriginalReads=true" \
	--export "--preset-file $PRESET" \
	$FASTQ \
	$EXPORT

### Move output files
mv $EXPORT\.clna $OUT/clna/
mv $EXPORT\.vdjca $OUT/vdjca/
mv $EXPORT\.clonotypes.TRB.txt $OUT/export_clones
