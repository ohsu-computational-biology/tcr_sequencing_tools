#!/usr/bin/sh

#SBATCH --partition		exacloud
#SBATCH --nodes			1
#SBATCH --ntasks		1
#SBATCH --ntasks-per-core	1
#SBATCH --cpus-per-task		1
#SBATCH --mem-per-cpu		16000
#SBATCH --output		homeo-%j.out
#SBATCH --error			homeo-%j.err

MYBIN=$tool/analysis/homeo
IN=$data/normalization/normalized_clones
OUT=$data/freqGroups
META=$data/QC/meta/LIB170920LC_meta.txt

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_JOB_NODELIST: " $SLURM_JOB_NODELIST
echo "SLURM_CPUS_ON_NODE: " $SLURM_CPUS_ON_NODE
echo "Current file: " $CURRFILE

mkdir -p $OUT

$MYBIN -i $IN -o $OUT -m $META -l TRUE

