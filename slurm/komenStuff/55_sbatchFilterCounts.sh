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
#SBATCH --output             filterClones_%A_%a.out        # Standard output
#SBATCH --error              filterClones_%A_%a.err        # Standard error
#SBATCH --array              0-5

### Dirs to run
DIRS=(LIB170728LC LIB170921LC LIB171213LC LIB190603LC LIB190701LC LIB191025LC)

### Arguments
CUTOFF=5
OLD=F

### Tool
MYBIN=$tool/40_postProcess/filterClones.R

## Get current
DIR=${DIRS[$SLURM_ARRAY_TASK_ID]}
data=$dha/$DIR

## Set vars
IN=$data/normalization/normalized_clones/
OUT=$data/normalization/filtered_clones_komen_$CUTOFF
QC=$data/QC/komen/

### Make directories
mkdir -p $QC
mkdir -p $OUT

cmd="/usr/bin/Rscript $MYBIN --inputDir $IN --outputDir $OUT --cutOff $CUTOFF --qcDir $QC --old $OLD" 

echo $cmd
eval $cmd
