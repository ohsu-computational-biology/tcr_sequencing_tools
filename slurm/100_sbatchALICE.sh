#!/bin/bash

#SBATCH --partition          exacloud                # partition (queue)
#SBATCH --nodes              1                       # number of nodes
#SBATCH --ntasks             1                       # number of "tasks" to be allocated for the job
#SBATCH --ntasks-per-core    1                       # Max number of "tasks" per core.
#SBATCH --cpus-per-task      4                       # Set if you know a task requires multiple processors
#SBATCH --mem-per-cpu        8000                    # Memory required per allocated CPU (mutually exclusive with mem)
##SBATCH --mem                16000                  # memory pool for each node
#SBATCH --time               0-24:00                 # time (D-HH:MM)
#SBATCH --output             alice_%A_%a.out        # Standard output
#SBATCH --error              alice_%A_%a.err        # Standard error
#SBATCH --array              1-7                    # sets number of jobs in array


### SET I/O VARIABLES
IN="$data/ALICE/data"
OUT="$data/ALICE/output"
LOAD="$lustre/myApps/ALICE"
CPU=4
N=1000000
MYBIN="$tool/70_advAnalysis/runAlice.R"          

### Record slurm info
date
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
printf "\n\n"

### Get file/directory to run on
#TODO=$data/tools/todo/alice.txt
#CURRDIR=`awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}' $TODO`
CURRDIR=`ls -v $IN | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

### Check inputs
printf "Main input directory: %s\n" $IN
printf "\tCurrent directory within main: %s\n" $CURRDIR
printf "Main output directory: %s\n" $OUT
printf "Load directory: %s\n" $LOAD
printf "\n\n"

### Execute
mkdir -p $OUT
cmd="/usr/bin/Rscript $MYBIN --inputDir $IN/$CURRDIR --outDir $OUT/$CURRDIR --loadDir $LOAD --cpu $CPU --nSeq $N"

echo $cmd
eval $cmd

