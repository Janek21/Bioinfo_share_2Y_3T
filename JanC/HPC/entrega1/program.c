#include <mpi.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv); // We initialize the MPI environment

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size); // Get the number of processes

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); // Get the current process

    char broadcast_message[100] = "Helloo :)";
    // Broadcast: The root process sends the message to all processes
    MPI_Bcast(broadcast_message, 100, MPI_CHAR, 0, MPI_COMM_WORLD);

    // Display the message received from the broadcast indicating the process
    printf("Process %d received broadcast message: %s\n", world_rank, broadcast_message);

    // We include the process ID
    char response_message[120];
    sprintf(response_message, "Process %d processed the message.", world_rank);

    // If we are not the root process ( root = world_rank = 0), send the modified message back to the root
    if (world_rank != 0) {
        MPI_Send(response_message, strlen(response_message) + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    } else {
        // The root process receives all messages 
        for (int i = 1; i < world_size; i++) {
            char recv_buf[120];
            MPI_Recv(recv_buf, 120, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Root process received: %s\n", recv_buf);
        }
    }

    MPI_Finalize();
    return 0;
}

