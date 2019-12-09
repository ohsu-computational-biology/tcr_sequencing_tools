#!/usr/bin/sh

#SBATCH --partition		exacloud
#SBATCH --nodes			1
#SBATCH --ntasks		1
#SBATCH --ntasks-per-core	1
#SBATCH --cpus-per-task		1
#SBATCH --mem-per-cpu		16000
#SBATCH --output		topCloneFreq-%j.out
#SBATCH --error			topCloneFreq-%j.err

MYBIN=$tool/60_analysis/topCloneFreq.R
IN=$data/normalization/normalized_clones/
OUT=$data/QC/std
FREQ="5,10,15,20,50,100,200"

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_JOB_NODELIST: " $SLURM_JOB_NODELIST
echo "SLURM_CPUS_ON_NODE: " $SLURM_CPUS_ON_NODE

mkdir -p $OUT

cmd="$MYBIN -i $IN -f $FREQ -o $OUT"

echo $cmd
eval $cmd

