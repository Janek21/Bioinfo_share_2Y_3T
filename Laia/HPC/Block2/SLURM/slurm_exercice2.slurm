#!/bin/bash
#SBATCH --job-name=array-job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=nodo.q
#SBATCH --array=0-4

# GET THE CORRESPONDENT FILE
Taskid=$SLURM_ARRAY_TASK_ID
file=$(ls names_*.txt | sed -n "$(($Taskid + 1))p")  #Get the correspondent file

# TRANSFORM THE NAME AND THE CONTENT OF THE FILE TO UPPERCASE
newname=$(basename "$file" .txt)
uppername=$(echo "$newname" | tr '[:lower:]' '[:upper:]')
tr '[:lower:]' '[:upper:]' < "$file" > "${uppername}.txt"

# SORT THE CONTENT OF THE FILE
sort "${uppername}.txt" -o "${uppername}.txt"

# MAKE A SUMMARY
./summary.sh "${uppername}.txt" > "summary_$Taskid.txt"

# MERGE SUMMARIES
cat "summary_$Taskid.txt" >> summary.txt

# REMOVE UNNECESSARY FILES
rm "summary_$Taskid.txt"q

