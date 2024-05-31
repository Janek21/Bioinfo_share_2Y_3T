#!/bin/bash
#SBATCH --job-name=convert_to_upper
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=nodo.q


for file in names_*.txt; do
  newname=$(basename "$file" .txt)
  uppername=$(echo "$newname" | tr '[:lower:]' '[:upper:]')
  
  tr '[:lower:]' '[:upper:]' < "$file" > "${uppername}.txt"
done

