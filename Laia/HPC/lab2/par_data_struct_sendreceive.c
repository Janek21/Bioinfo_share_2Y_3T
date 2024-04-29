#include <stdio.h>
#include <unistd.h>
#include "mpi.h"

/* This code uses illustrates MPI_Sendrecv */

/* This example handles a 12 x 12 mesh processed on P=4 processes. */
/* NOTE: The program WILL ABORT if NOT RUN WITH EXACTLY P PROCESSES. */

/* HOMEWORK:
     Fill in the missing parameters in Statement Labelled S12 */


#define maxn 12
#define P 4

int main( argc, argv )
int argc;
char **argv;
{
    int rank, value, size, errcnt, toterr, i, j;
    int next_nbr, prev_nbr;
    MPI_Status status;
    //double x[maxn][maxn];
    double xlocal[(maxn/P)+2][maxn];

    MPI_Init( &argc, &argv );

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    if (size != P) {
       printf("This program only works when run with %d processes. Aborting.\n", P);
       usleep(10);
       MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    /* xlocal[0][] refers to the initial row. It's used to store ghost points
                   used to keep data coming from the previous process */

    /* xlocal[maxn/size+1][] refers to the final row: it's used to store ghost
                   points used to keep data coming from the next process */

    /* Fill the data as specified */
    /* 1) Fill its local segment with the rank of the process */
    for (i=1; i<=maxn/size; i++) 
	for (j=0; j<maxn; j++) 
	    xlocal[i][j] = rank;
    /* 2) Fill ghost points with -1 */
    for (j=0; j<maxn; j++) {
	xlocal[0][j] = -1;
	xlocal[maxn/size+1][j] = -1;
    }

    /* We use xlocal[i] for &xlocal[i][0] to provide the initial address
       of the i-th row, since both expressions give the same memory address */

    /* Note that we use MPI_PROC_NULL to remove the if statements that
       would be needed without MPI_PROC_NULL */
    next_nbr = rank + 1;
    if (next_nbr >= size) next_nbr = MPI_PROC_NULL;
    prev_nbr = rank - 1;
    if (prev_nbr < 0) prev_nbr = MPI_PROC_NULL;

    /* Send last local row to next process (to process with identifier rank+1)
       unless I'm the last one (process with identifier size-1)
                               &&
       Receive from previous unless I'm the first process (rank==0)
       and store row in initial row (ghost area) */
    MPI_Sendrecv( xlocal[maxn/size], maxn, MPI_DOUBLE, next_nbr, 0,
		  xlocal[0], maxn, MPI_DOUBLE, prev_nbr, 0, 
		  MPI_COMM_WORLD, &status );

    /* Send 1st local row to previous process unless I'm the 1st process 
                               &&
       Receive from next process and store row in last row in xlocal
       (ghost area) unless I'm the last process */
    MPI_Sendrecv( xlocal[1], maxn, MPI_DOUBLE, prev_nbr, 1, 
		  xlocal[maxn/size+1], maxn, MPI_DOUBLE, next_nbr, 1, 
		  MPI_COMM_WORLD, &status );   /* Statement S12 */

    /* Check that we have the correct results */
    errcnt = 0;
    for (i=1; i<=maxn/size; i++) 
	for (j=0; j<maxn; j++) 
	    if (xlocal[i][j] != rank) errcnt++;
    for (j=0; j<maxn; j++) {
	if (xlocal[0][j] != rank - 1) errcnt++;
	if (rank < size-1 && xlocal[maxn/size+1][j] != rank + 1) errcnt++;
    }

    /* Retrieve in process 0 the total number of errors */
    MPI_Reduce( &errcnt, &toterr, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
	if (toterr)
	    printf( "ERROR: found %d errors\n", toterr );
	else
	    printf( "No errors\n" );
    }

    MPI_Finalize( );
    return 0;
}
