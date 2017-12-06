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
#SBATCH --output             subsetClones_%j.out        # Standard output
#SBATCH --error              subsetClones_%j.err        # Standard error


### SET I/O VARIABLES

IN1=$data/normalization/collapsed_clones/             # Directory containing all input files. Should be one job per file
IN2=$data/freqGroups/collapse_groupData/newNorm/LIB170920LC_full_clones.txt
GRP=
OUT=$data/vdjtools/           # Directory where output files should be written
QC=$data/QC/
MYBIN=$tool/misc/subsetByGroup.R          # Path to shell script or command-line executable that will be used

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

GRPS=('Hyperexpanded,Large' 'Hyperexpanded,Large,Medium' 'Medium' 'Rare,Small' 'Rare' 'Hyperexpanded,Large,Medium,Small,Rare')
DIRS=('hyperLarge' 'hyperLargeMedium' 'medium' 'rareSmall' 'rare' 'all')

for i in {0..6}; do

    ## Get variables
    CURRGRP=${GRPS[$i]}
    CURRDIR=${DIRS[$i]}

    ## Update
    echo "Current group designation: " $CURRGRP
    echo "Current output directory: " $OUT/$CURRDIR/

    cmd="/usr/bin/Rscript $MYBIN -i $IN1 -f $IN2 -d $CURRGRP -o $OUT/$CURRDIR/ -l F" 

    echo $cmd
    eval $cmd
done
