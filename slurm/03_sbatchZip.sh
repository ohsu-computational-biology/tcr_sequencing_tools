#!/bin/bash

### Unzip all of the fastq files. Need to first determine the number of different tens place file names there are. There will be that many array jobs.

#SBATCH --partition          exacloud                # partition (queue)
#SBATCH --nodes              1                       # number of nodes
#SBATCH --ntasks             1                       # number of "tasks" to be allocated for the job
#SBATCH --ntasks-per-core    1                       # Max number of "tasks" per core.
#SBATCH --cpus-per-task      1                       # Set if you know a task requires multiple processors
#SBATCH --mem                16000                  # memory pool for each node
#SBATCH --time               0-24:00                 # time (D-HH:MM)
#SBATCH --output             zip_%A_%a.out        # Standard output
#SBATCH --error              zip_%A_%a.err        # Standard error
#SBATCH --array              0-5                    # sets number of jobs in array


### SET I/O VARIABLES

BATCH=LIB170920LC
IN=$dha/$BATCH/fastqs_from_core/decontam/ 	# Directory containing all input files. Should be one job per file
MYBIN=$tool/misc/zip.sh                    # Path to shell script or command-line executable that will be used

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

### Get a template file
TEMP=`ls -v $IN | head -1`
BASE="${TEMP%%S[0-9]*}"
TENS=$SLURM_ARRAY_TASK_ID

### Get suffix
SUFFIX="${TEMP#$BATCH\_S[0-9]*}"
FASTQ=`echo $SUFFIX | grep _R2`

if [ -v "$FASTQ" ]; then
    SUFFIX2=""
else
    SUFFIX2=`echo $SUFFIX | sed 's/R1/R2/'`
fi


### Print checks
echo "Example file: " $TEMP
echo "File base: " $BASE
echo "Suffix 1: " $SUFFIX
echo "Suffix 2: " $SUFFIX2
printf "\n\n"

### Execute

for i in {0..9}; do

    ## Get file names
    if [ $TENS == 0 ]; then
        FILE1=$IN/$BASE\S$i$SUFFIX
        FILE2=$IN/$BASE\S$i$SUFFIX2
    else
        FILE1=$IN/$BASE\S$TENS$i$SUFFIX
        FILE2=$IN/$BASE\S$TENS$i$SUFFIX2
    fi

    ## Check file names
    echo "Files to run: " $FILE1 $FILE2

    ## Prepare command
    cmd1="sh $MYBIN $FILE1"
    cmd2="sh $MYBIN $FILE2"

    ## Echo/Evaluate command 1
    echo $cmd1
    eval $cmd1

    ## Echo/Evaluate command 2
    if [ ! -v $SUFFIX2 ]; then
        echo $cmd2
        eval $cmd2
    fi

    printf "\n\n"

done
