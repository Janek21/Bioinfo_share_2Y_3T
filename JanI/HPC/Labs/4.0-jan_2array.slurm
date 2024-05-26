#!/bin/bash
#SBATCH --job-name=uppercase_sumup-job
#SBATCH --nodes=5
#SBATCH --ntasks=5
#SBATCH --ntasks-per-node=1
#SBATCH --partition=nodo.q
#SBATCH --array=0-4


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

rm summary_$Tid.txt #unclog the directory

# print initial date and time
echo "End UPPERCASE JOB at $(date)"

#If we were to use another sbatch command we would use it in the following manner (without the #)
#sbatch submit_example.<slurm or sh>
#for example:
#sbatch summary.sh names_*.txt
#this would get the summary results in a file named slurm-<jobID>.out
#and to get the job id we are to use $SLURM_JOB_ID as a variable
#to create a dependending execution the command would be
#sbatch --dependency=<ID of the one the next task depends on> submit_example.<slurm or sh>
