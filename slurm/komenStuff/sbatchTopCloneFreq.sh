#!/usr/bin/sh

#SBATCH --partition		exacloud
#SBATCH --nodes			1
#SBATCH --ntasks		1
#SBATCH --ntasks-per-core	1
#SBATCH --cpus-per-task		1
#SBATCH --mem-per-cpu		16000
#SBATCH --output		topCloneFreq-%j.out
#SBATCH --error			topCloneFreq-%j.err

DIRS=(LIB170728LC LIB170921LC LIB171213LC LIB190603LC LIB190701LC LIB191025LC)

GROUP='"5,6-25,26-50,51-100"'
REST=T

for DIR in "${DIRS[@]}"; do

	data=$dha/$DIR
	
	MYBIN=$tool/60_analysis/topCloneFreq.R
	IN=$data/normalization/normalized_clones/
	OUT=$data/QC/std
	
	mkdir -p $OUT
	
	printf "Batch: %s\nGroups: %s\nOutput: %s\n" $DIR $GROUP $OUT
	
	cmd="$MYBIN -i $IN -g $GROUP -r $REST -o $OUT"
	
	echo $cmd
	eval $cmd
done

