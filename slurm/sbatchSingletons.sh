#!/usr/bin/sh

#SBATCH --partition		exacloud
#SBATCH --nodes			1
#SBATCH --ntasks		1
#SBATCH --ntasks-per-core	1
#SBATCH --cpus-per-task		1
#SBATCH --mem-per-cpu		16000
#SBATCH --output		singleton-%j.out
#SBATCH --error			singleton-%j.err

DIRS=(LIB170728LC LIB170921LC LIB171213LC LIB190603LC LIB190701LC LIB191025LC)

for DIR in "${DIRS[@]}"; do

	data=$dha/$DIR
	
	MYBIN=$tool/60_analysis/countSingletons.R
	IN=$data/normalization/normalized_clones/
	OUT=$data/QC/std
	
	mkdir -p $OUT
	
	cmd="$MYBIN -i $IN -o $OUT"
	
	echo $cmd
	eval $cmd
done

