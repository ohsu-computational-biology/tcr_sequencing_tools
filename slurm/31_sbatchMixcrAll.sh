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
#SBATCH --mem                32000                  # memory pool for each node
#SBATCH --time               0-24:00                 # time (D-HH:MM)
#SBATCH --output             complete_mixcr_%A_%a.out        # Standard output
#SBATCH --error              complete_mixcr_%A_%a.err        # Standard error
#SBATCH --array              1-1                    # sets number of jobs in array


### SET I/O VARIABLES

IN=$data/mixcr/despiked_fastqs             # Directory containing all input files. Should be one job per file
OUT=$data/mixcr/           # Directory where output files should be written
REPORT=$data/mixcr/reports/
INDEX=$data/mixcr/indexes
MYBIN=$MIXCR/mixcr.jar          # Path to shell script or command-line executable that will be used
PRESET=$tool/mixcr/assemble_preset.txt
SPECIES="mmu"
TEST=False

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

#################
### VARIABLES ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################

### todo file
TODO=$data/tools/todo/align.txt
CURRFILE=`awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}' $TODO`
BASE="${CURRFILE%%S[0-9]*}"
SNUM="${CURRFILE%%.assembled.removed.fastq}"; SNUM="${SNUM##*S}"

### Issign I/O files/paths
FASTQ=$IN/$CURRFILE
ALIGNREPORT=S$SNUM\_align_report.txt
ALIGNMENT=$BASE\S$SNUM\_alignment.vdjca

ASSEMBLEREPORT=S$SNUM\_assemble_report.txt
ASSEMBLY=$BASE\S$SNUM\_clones.clns
INDEXFILE=$INDEX/S$SNUM\_index

EXPORT=$BASE/S$SNUM\_clones.txt

#############
### ALIGN ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############

cmd="/usr/bin/java -Xmx15g -jar $MYBIN align -f -c TRB --species $SPECIES --save-description --save-reads -v --report $REPORT/align/$ALIGNREPORT $FASTQ $OUT/alignments/$ALIGNMENT"

#cmd="/usr/bin/java -Xmx15g -jar $MYBIN align -f -c TRB --species mmu --save-description --save-reads -v --report $REPORTFILE $FASTQ $OUTPUT"

echo "ALIGN"
echo $cmd
printf "\n\n"
if [ $TEST != True ]; then eval $cmd; else echo "Testing"; fi

################
### ASSEMBLE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################

cmd="/usr/bin/java -Xmx15g -jar $MYBIN assemble -f --index $INDEXFILE --report $REPORT/assemble/$ASSEMBLEREPORT --threads 4 $OUT/alignments/$ALIGNMENT $OUT/assemblies/$ASSEMBLY"

echo "ASSEMBLE"
echo $cmd
printf "\n\n"
if [ $TEST != True ]; then eval $cmd; else echo "Testing"; fi

##############
### EXPORT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############

cmd="/usr/bin/java -Xmx15g -jar $MYBIN exportClones -f --preset-file $PRESET -o -t -readIds $INDEXFILE -cloneId $INDEXFILE $OUT/assemblies/$ASSEMBLY $OUT/export_clones/$EXPORT"

echo "EXPORT"
echo $cmd
printf "\n\n"
if [ $TEST != True ]; then eval $cmd; else echo "Testing"; fi


