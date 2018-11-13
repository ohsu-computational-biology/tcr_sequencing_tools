#!/usr/bin/sh

#SBATCH --partition		exacloud
#SBATCH --nodes			1
#SBATCH --ntasks		1
#SBATCH --ntasks-per-core	1
#SBATCH --cpus-per-task		1
#SBATCH --mem-per-cpu		16000
#SBATCH --output		clonalDiv-%j.out
#SBATCH --error			clonalDiv-%j.err

MYBIN=$tool/60_analysis/clonalDivisionSummary.R
IN=$data/normalization/normalized_clones/
OUT=$data/QC/cloneDiv/
META=$data/QC/meta/meta.txt

### Commands
### 1. cloneDir_v
### 2. metadata_v
### 3. outDir_v
### 4. toWrite_v
### 5. old_v
### 6. tissue_v (optional)
### 7. type_v (optional, but requires tissue_v)

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_JOB_NODELIST: " $SLURM_JOB_NODELIST
echo "SLURM_CPUS_ON_NODE: " $SLURM_CPUS_ON_NODE
echo "Current file: " $CURRFILE

mkdir -p $OUT

TOWRITE=TRUE
OLD=FALSE

cmd="$MYBIN \
	--cloneDir $IN \
	--metadata $META \
	--outDir $OUT \
	--write $TOWRITE \
	--old $OLD" 
	#--tissue $TISSUE \
	#--type $TYPE"

echo $cmd
eval $cmd

