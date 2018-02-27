#!/bin/bash


#SBATCH --partition          exacloud                # partition (queue)
#SBATCH --nodes              1                       # number of nodes
#SBATCH --ntasks             1                       # number of "tasks" to be allocated for the job
#SBATCH --ntasks-per-core    1                       # Max number of "tasks" per core.
#SBATCH --cpus-per-task      1                       # Set if you know a task requires multiple processors
##SBATCH --mem-per-cpu        8000                    # Memory required per allocated CPU (mutually exclusive with mem)
#SBATCH --mem                32000                  # memory pool for each node
#SBATCH --time               0-24:00                 # time (D-HH:MM)
#SBATCH --output             complete_mixcr_RNA_%A_%a.out        # Standard output
#SBATCH --error              complete_mixcr_RNA_%A_%a.err        # Standard error
#SBATCH --array              1-1                    # sets number of jobs in array


### SET I/O VARIABLES

IN=$data/fastqs_from_core/fastqs/             # Directory containing all input files. Should be one job per file
OUT=$data/mixcr/alignments        # Directory where output files should be written
REPORT=$data/mixcr/reports/align
INDEX=$data/mixcr/indexes
MIXCR=/home/exacloud/lustre1/CompBio/users/hortowe/tcr_sequencing_tools/mixcr/mixcr-2.1.9
MYBIN=$MIXCR/mixcr.jar          # Path to shell script or command-line executable that will be used
PRESET=$tool/mixcr/assemble_preset.txt
SPECIES="mmu"                   # mmu for mouse or hsa for human
TEST=False			# set to True to print out commands instead of running them

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

### Todo file
TODO=$data/tools/todo/align.txt
CURRFILE=`awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}' $TODO`
BASE="${CURRFILE%%S[0-9]*}"
SNUM="${CURRFILE%%.fastq}"; SNUM="${SNUM##*S}"

### Issign I/O files/paths

FASTQ=$IN/$CURRFILE

ALIGNREPORT=S$SNUM\_align_report.txt
ALIGNMENT=$BASE\S$SNUM\_alignment.vdjca

PARTIAL1=$BASE\S$SNUM\_partial1.vdjca
PARTIAL1REPORT=S$SNUM\_partial1_report.txt

PARTIAL2=$BASE\S$SNUM\_partial2.vdjca
PARTIAL2REPORT=S$SNUM\_partial2_report.txt

EXTENDED=$BASE\S$SNUM\_extended.vdjca

ASSEMBLEREPORT=$REPORT/S$SNUM\_assemble_report.txt
ASSEMBLY=$BASE\S$SNUM\_clones.clns
INDEXFILE=$INDEX/S$SNUM\_index.txt

EXPORT=$BASE/S$SNUM\_clones.txt

#############
### ALIGN ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############

cmd="/usr/bin/java -Xmx15g -jar $MYBIN align -f --species $SPECIES -OallowPartialAlignments=true --save-description --save-reads -v --report $REPORT/$ALIGNREPORT $FASTQ $OUT/$ALIGNMENT"

echo "ALIGN"
echo $cmd
printf "\n\n"
if [ $TEST != True ]; then eval $cmd; else echo "Testing"; fi

################
### PARTIAL1 ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################

### New directories
IN=$data/mixcr/alignments
OUT=$data/mixcr/partial1
REPORT=$data/mixcr/reports/partial1

cmd="/usr/bin/java -Xmx15g -jar $MYBIN assemblePartial -f -r $REPORT/$PARTIAL1REPORT $IN/$ALIGNMENT $OUT/$PARTIAL1"

echo "PARTIAL1"
echo $cmd
printf "\n\n"
if [ $TEST != True ]; then eval $cmd; else echo "Testing"; fi

################
### PARTIAL2 ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################

### New directories
IN=$data/mixcr/partial1
OUT=$data/mixcr/partial2
REPORT=$data/mixcr/reports/partial2

cmd="/usr/bin/java -Xmx15g -jar $MYBIN assemblePartial -f -r $REPORT/$PARTIAL2REPORT $IN/$PARTIAL1 $OUT/$PARTIAL2"

echo "PARTIAL2"
echo $cmd
printf "\n\n"
if [ $TEST != True ]; then eval $cmd; else echo "Testing"; fi

##############
### EXTEND ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############

### New directories
IN=$data/mixcr/partial2
OUT=$data/mixcr/extended

cmd="/usr/bin/java -Xmx15g -jar $MYBIN extendAlignments -f $IN/$PARTIAL2 $OUT/$EXTENDED"

echo "EXTEND"
echo $cmd
printf "\n\n"
if [ $TEST != True ]; then eval $cmd; else echo "Testing"; fi

################
### ASSEMBLE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################

### New directories
IN=$data/mixcr/extended
OUT=$data/mixcr/assemblies
REPORT=$data/mixcr/reports/assemble

cmd="/usr/bin/java -Xmx15g -jar $MYBIN assemble -f --index $INDEXFILE --report $REPORT/$ASSEMBLEREPORT --threads 4 $IN/$EXTENDED $OUT/$ASSEMBLY"

echo "ASSEMBLE"
echo $cmd
printf "\n\n"
if [ $TEST != True ]; then eval $cmd; else echo "Testing"; fi


##############
### EXPORT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############

### New directories
IN=$data/mixcr/assemblies
OUT=$data/mixcr/export_clones
REPORT=$data/mixcr/reports/assemble

cmd="/usr/bin/java -Xmx15g -jar $MYBIN exportClones -f --preset-file $PRESET -o -t -readIds $INDEXFILE -cloneId $INDEXFILE $IN/ASSEMBLY $OUT/$EXPORT"

echo "EXPORT"
echo $cmd
printf "\n\n"
if [ $TEST != True ]; then eval $cmd; else echo "Testing"; fi


