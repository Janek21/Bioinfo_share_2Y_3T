#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv) {
    int rank, size;	// define the rank and the size as integers
    MPI_Init(&argc, &argv);	// initialize the MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);	// get the rank
    MPI_Comm_size(MPI_COMM_WORLD, &size);	// get the size

    int N = 10;	// number up to which we want to sum.
    int local_sum = 0; // local sum (for each process)
    int total_sum = 0; // global sum (sum for all the processes)

    // each process computes a part of the sum
    for (int i = rank + 1; i <= N; i += size) {
        local_sum += i;
    }

    // each process sends its local_sum to the root process (rank 0)
    if (rank != 0) {
        MPI_Send(&local_sum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    } else {
        total_sum = local_sum; // for root process, total_sum is initialized with its local_sum
        for (int p = 1; p < size; p++) {
            MPI_Recv(&local_sum, 1, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_sum += local_sum; // add the local_sum received from other processes
        }
    }

    // the process with rank = 0 will print out the result
    if (rank == 0) {
        printf("The sum of the first %d numbers is: %d\n", N, total_sum);
    }

    MPI_Finalize();	// finalize the execution
    return 0;
}
	
