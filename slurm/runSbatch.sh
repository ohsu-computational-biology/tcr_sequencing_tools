#!/bin/sh

### Set variables
SCRIPT=$1        # sbatch script to be submitted
FILES=$2         # maximum number of files

### Going to run for every 10 files
TENS=$(expr $FILES / 10)

for i in $(seq 0 $TENS); do

    counter=0
    ## Replace the array indexes
    if [ $i == 0 ]; then
	sed -E -i 's/(\#SBATCH --array)(.*$)/\1 1-9/' $SCRIPT
    else
	first=$i\0
	second=$i\9
	sed -E -i "s/(\#SBATCH --array)(.*$)/\1 $first-$second/" $SCRIPT
    fi

    echo "Currently on: " $i

    ## Run sbatch and capture
    curr_sub=$(sbatch $SCRIPT)
    curr_proc=${curr_sub##*\ }

    echo "Process name: " $curr_proc

    ## Wait till this one is done before going on to the next one
    #echo "Before while"
    while true; do
	curr_status=`eval squeue -u hortowe | cut -d ' ' -f 9 | grep $curr_proc`
	if [ -z "$curr_status" ]; then
	    echo "Done with: " $i
	    break
	elif [ ${#curr_status} == 1 ]; then
	    counter=$((counter + 1))
	    sleep 10
	    if [ counter >= 30 ]; then
		echo "Failed runs: " $curr_status
		scancel $curr_proc
	    fi
	else
	    :
	fi
    done
done
