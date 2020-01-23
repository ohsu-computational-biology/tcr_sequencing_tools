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
#SBATCH --output             divAnalysis_%A_%a.out        # Standard output
#SBATCH --error              divAnalysis_%A_%a.err        # Standard error
#SBATCH --array              0-5

### Dirs to run
DIRS=(LIB170728LC LIB170921LC LIB171213LC LIB190603LC LIB190701LC LIB191025LC)

### Cut-offs to run
CUTS=(1 5 10)

### SET I/O VARIABLES
DIR=${DIRS[$SLURM_ARRAY_TASK_ID]}
data=$dha/$DIR
BASE=$data/normalization/filtered_clones_komen_
OUT=$data/QC/komen/
MYBIN=$tool/60_analysis/diversityAnalysis.R
OLD=FALSE

mkdir -p $OUT

for CUT in ${CUTS[@]}; do

	IN=$BASE$CUT

	cmd="/usr/bin/Rscript $MYBIN $IN $OUT $OLD" 

	echo $cmd
	eval $cmd
	
	mv $OUT/$DIR\.txt $OUT/$DIR\_$CUT\.txt
done

