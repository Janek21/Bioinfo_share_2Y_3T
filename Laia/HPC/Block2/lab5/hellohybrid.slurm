#!/bin/bash

# A job submission script for running a hybrid MPI/OpenMP job
# Regard that not all the options might be needed.
# Just use the ones that are needed to get the right combination 
# of nodes, processes and threads

#SBATCH --job-name=hellohybrid
#SBATCH --output=hellohybrid_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=
#SBATCH --ntasks-per-socket=
#SBATCH --cpus-per-task=12
#SBATCH --partition=nodo.q

# Load the default OpenMPI module.

module load openmpi/4.1.4

# Set OMP_NUM_THREADS to the number of CPUs per task we asked for

export OMP_NUM_THREADS=12

# Run the process with mpirun. Note that the -n option is not required
# in this case; mpirun will automatically determine how many processes
# to run from the Slurm settings.

mpirun hellohybrid

