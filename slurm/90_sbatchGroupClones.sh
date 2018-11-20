#!/usr/bin/sh

#SBATCH --partition		exacloud
#SBATCH --nodes			1
#SBATCH --ntasks		1
#SBATCH --ntasks-per-core	1
#SBATCH --cpus-per-task		1
#SBATCH --mem-per-cpu		16000
#SBATCH --output		groupClones-%j.out
#SBATCH --error			groupClones-%j.err

MYBIN=$tool/60_analysis/groupClones.R
IN=$data/normalization/normalized_clones/
OUT=$data/freqGroups/
META=$data/QC/meta/LIB170920LC_meta.txt

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_JOB_NODELIST: " $SLURM_JOB_NODELIST
echo "SLURM_CPUS_ON_NODE: " $SLURM_CPUS_ON_NODE
echo "Current file: " $CURRFILE

mkdir -p $OUT

### Additional Arguments
COLUMNS="'nb.clone.count,nb.clone.fraction,nSeqCDR3,aaSeqCDR3,V segments,J segments'"                         # New norm, new mixcr column names

DIVISIONS="'Rare = 0.00001,Small = 0.0001,Medium = 0.001,Large = 0.01,Hyperexpanded = 1'"

LOG=TRUE

cmd="$MYBIN -i $IN -o $OUT -c $COLUMNS -m $META -d $DIVISIONS -l $LOG" 

echo $cmd
eval $cmd

