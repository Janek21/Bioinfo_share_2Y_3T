#!/bin/bash
#SBATCH --job-name=uppercase_sumup-job
#SBATCH --nodes=5
#SBATCH --ntasks=5
#SBATCH --ntasks-per-node=1
#SBATCH --partition=nodo.q
#SBATCH --array=0-4
# posar aixo enlloc del 4->$(ls names_*.txt | wc -l) #no es pot


# print initial date and time
echo "Start UPPERCASE JOB  at $(date)"
echo "-------------------------"
# print name of host
hostname
echo "-------------------------"

Tid=$SLURM_ARRAY_TASK_ID

f=names_$Tid.txt  #array 0-4, names_[0 to 4].txt, so each task will do a file

uppername=$(echo ${f::-4}| tr [:lower:] [:upper:]) #transform each name to upper case without the .txt

cat $f | tr [:lower:] [:upper:] > "$uppername.txt" # sort the contents of each file and turn them uppercase, place them in NAME_*.txt

# print initial date and time
echo "End UPPERCASE JOB at $(date)"
