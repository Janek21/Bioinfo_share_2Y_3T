#!/bin/bash
#SBATCH --job-name=array-job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=nodo.q
#SBATCH --array=0-4
#trobar com canviar el 4 pel numero de files

# print initial date and time
echo "Start UPPERCASE JOB  at $(date)"
echo "-------------------------"
# print name of host
hostname
echo "-------------------------"

Tid=$SLURM_ARRAY_TASK_ID

f=names_$Tid.txt  #array 0-4, names_[0 to 4].txt, so each task will do a file

uppername=$(echo ${f::-4}| tr [:lower:] [:upper:]) #transform each name to upper case without the .txt

sort $f | tr [:lower:] [:upper:] > "$uppername.txt" # sort the contents of each file and turn them uppercase, place them in NAME_*.txt

#use %1 in --array to execute 1 task a time
#sleep $Tid

./summary.sh "$uppername.txt" > summary_$Tid.txt #create the summary for each names file

cat summary_$Tid.txt >> summary.txt

rm summary_$Tid.txt

# print initial date and time
echo "End UPPERCASE JOB at $(date)"
