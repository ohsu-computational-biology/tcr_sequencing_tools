#!/bin/bash

### This script provides a template containg many of the commonly-used resource management SBATCH commands.
### You can submit this script as-is using `sbatch sbatchTemplate.sh` and it will output the slurm info into each log file.
### Use the following commands to combine the 10 (default) log files into a space-delimited file for testing purposes.
### From directory of output files (directory where sbatch is run from):

#SBATCH --partition          exacloud                # partition (queue)
#SBATCH --nodes              1                       # number of nodes
#SBATCH --ntasks             1                       # number of "tasks" to be allocated for the job
#SBATCH --ntasks-per-core    1                       # Max number of "tasks" per core.
#SBATCH --cpus-per-task      1                       # Set if you know a task requires multiple processors
##SBATCH --mem-per-cpu        8000                    # Memory required per allocated CPU (mutually exclusive with mem)
#SBATCH --mem                16000                  # memory pool for each node
#SBATCH --time               0-24:00                 # time (D-HH:MM)
#SBATCH --output             normalize_%A_%a.out        # Standard output
#SBATCH --error              normalize_%A_%a.err        # Standard error
#SBATCH --array              1-1                    # sets number of jobs in array


### SET I/O VARIABLES

IN=$data/normalization             # Directory containing all input files. Should be one job per file
REF=$tool/40_postProcess/nb_counts/
MYBIN=$tool/40_postProcess/normalizeClones.R          # Path to shell script or command-line executable that will be used

### Record slurm info

date
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
echo "SLURM_JOB_NODELIST: " $SLURM_JOB_NODELIST
echo "SLURM_CPUS_ON_NODE: " $SLURM_CPUS_ON_NODE
echo "SLURM_CPUS_PER_TASK: " $SLURM_CPUS_PER_TASK
echo "SLURM_JOB_CPUS_PER_NODE: " $SLURM_JOB_CPUS_PER_NODE
echo "SLURM_MEM_PER_CPU: " $SLURM_MEM_PER_CPU
echo "SLURM_MEM_PER_NODE: " $SLURM_MEM_PER_NODE
echo "SLURM_NTASKS: " $SLURM_NTASKS
echo "SLURM_NTASKS_PER_CORE " $SLURM_NTASKS_PER_CORE
echo "SLURM_NTASKS_PER_NODE " $SLURM_NTASKS_PER_NODE
echo "SLURM_TASKS_PER_NODE " $SLURM_TASKS_PER_NODE
printf "\n\n"

### create array of file names in this location (input files)
### Extract file name information as well

TODO=$data/tools/todo/normalize.txt
CURRFILE=`awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}' $TODO`
#CURRFILE=`ls -v $IN/decontam | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
BASE="${CURRFILE%%S[0-9]*}"
SNUM="${CURRFILE%%_alignment_clones_exported_decontam.txt}"; SNUM="${SNUM##*S}"

echo "Clone File: " $CURRFILE
echo "Base name: " $BASE
echo "Sample number: " $SNUM
printf "\n\n"

### Set I/O files/paths
RAWCLONE=$IN/decontam/$CURRFILE
SPIKECOUNTS=$IN/counts/$BASE\S$SNUM\.assembled.spike.counts.25bp.txt
NORMDIR=$IN/normalized_clones/
ORIGSCALE=$IN/scaling_factor.txt
NBSCALE=$REF/nb.scaling.factors.txt

echo "Raw Count input: " $RAWCLONE
echo "Spike Count input: " $SPIKECOUNTS
echo "Normalized Count Output Directory: " $NORMDIR
echo "Original method scaling factors " $ORIGSCALE
echo "NB method scaling factors " $NBSCALE
printf "\n\n"

### Execute

cmd="/usr/bin/Rscript $MYBIN $RAWCLONE $SPIKECOUNTS $NORMDIR $ORIGSCALE $NBSCALE"

echo $cmd
eval $cmd

