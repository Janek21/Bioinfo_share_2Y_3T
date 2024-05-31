#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "mpi.h"

/* This example handles a 12 x 12 mesh processed on P=4 processes. */
/* NOTE: The program WILL ABORT if NOT RUN WITH EXACTLY P PROCESSES. */

#define maxn 12
#define P 4

int main( argc, argv )
int argc;
char **argv;
{
    int         rank, value, size, errcnt, toterr, i, j, itcnt;
    int         i_first, i_last;
    int         next_nbr, prev_nbr;
    int         numreq;
    MPI_Request r[4];
    MPI_Status  statuses[2];
    MPI_Status  status;
    double      diffnorm, globaldiffnorm;
    double      xlocal[(maxn/P)+2][maxn];
    double      xnew[(maxn/P)+1][maxn];     /* Auxiliary matrix used in the Jacobi method. */
    double      x[maxn][maxn];

    MPI_Init( &argc, &argv );

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    if (size != P) {
       printf("This program only works when run with %d processes. Aborting.\n", P);
       usleep(10);
       MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    /* xlocal[0][] refers to the row used to store ghostpoints
                   used to keep data coming from the previous process */

    /* xlocal[maxn/size+1][] refers to the row used to store ghostpoints
                   used to keep data coming from the next process */

    i_first = 1;
    i_last  = maxn/size;
    /* Note that first and last processes have one less row of interior
       points */
    if (rank == 0)        i_first++;
    if (rank == size - 1) i_last--;

    /* Fill the data as specified */
    for (i=1; i<=maxn/size; i++) 
	for (j=0; j<maxn; j++) 
	    xlocal[i][j] = rank;
    for (j=0; j<maxn; j++) {
	xlocal[i_first-1][j] = -1;
	xlocal[i_last+1][j] = -1;
    }


    /* Note that we use MPI_PROC_NULL to remove the if statements that
       would be needed without MPI_PROC_NULL */
    next_nbr = rank + 1;
    if (next_nbr >= size) next_nbr = MPI_PROC_NULL;
    prev_nbr = rank - 1;
    if (prev_nbr < 0) prev_nbr = MPI_PROC_NULL;

    itcnt = 0;
    numreq=2; /* Prepare to synchronize at once with two communication requests:
               - r[0]: sending last local row xlocal[maxn/size]
	       - r[1]: receiving ghost row and storing in xlocal[maxn/size+1] */
    do {
	/* Note the use of xlocal[i] for &xlocal[i][0] */

        /* Send last local row to next process (process with identifier rank+1)
           unless I'm the last one (process with identifier size-1). */
	MPI_Isend(xlocal[maxn/size],maxn,MPI_DOUBLE,next_nbr,0,MPI_COMM_WORLD,
                  &r[0] ); 
        /* Receive from previous unless I'm the first process (rank==0)
           and store row in initial row (ghost area) */
	MPI_Irecv(xlocal[0],maxn,MPI_DOUBLE,prev_nbr,0,MPI_COMM_WORLD,
                  &r[2] );

        /* Send 1st local row to previous process unless I'm the 1st process */
	MPI_Isend(xlocal[1],maxn,MPI_DOUBLE,prev_nbr,1,MPI_COMM_WORLD,
	          &r[3] );
        /* Receive from next process and store row in last row in xlocal
           (ghost area) unless I'm the last process */
	MPI_Irecv(xlocal[maxn/size+1],maxn,MPI_DOUBLE,next_nbr,1,MPI_COMM_WORLD,
	          &r[1] );
	
	/* Make sure we do NOT use row xlocal[0] before it is received */
	MPI_Wait( &r[2], MPI_STATUS_IGNORE ); /* Statement S16 (Fix the error) */
	/* Compute new values (but not on boundary) */
	itcnt ++;
	diffnorm = 0.0;
	for (i=i_first; i<i_last; i++) /* Note the < instead of <= */
	    for (j=1; j<maxn-1; j++) {
		xnew[i][j] = (xlocal[i][j+1] + xlocal[i][j-1] +
			      xlocal[i+1][j] + xlocal[i-1][j]) / 4.0;
		diffnorm += (xnew[i][j] - xlocal[i][j]) * 
		            (xnew[i][j] - xlocal[i][j]);
	    }

	/* Last iteration has been taken out of the loop */
	i=i_last; 
        /* Need to synchronize with pending communication requests envolving
           the last 2 rows. Synchronize at once with 2 communication requests:
               - r[0]: sending last local row xlocal[maxn/size]
	       - r[1]: receiving ghost row and storing in xlocal[maxn/size+1] */
	MPI_Waitall( numreq, r, MPI_STATUS_IGNORE ); /* Statement S17 */
	for (j=1; j<maxn-1; j++) {
		xnew[i][j] = (xlocal[i][j+1] + xlocal[i][j-1] +
			      xlocal[i+1][j] + xlocal[i-1][j]) / 4.0;
		diffnorm += (xnew[i][j] - xlocal[i][j]) * 
		            (xnew[i][j] - xlocal[i][j]);
	}

	/* Make sure we do not overwrite row xlocal[1] before it is sent */
	MPI_Wait( &r[3], MPI_STATUS_IGNORE ); /* Statement S18 */
	/* Only transfer the interior points */
	for (i=i_first; i<=i_last; i++) 
	    for (j=1; j<maxn-1; j++) 
		xlocal[i][j] = xnew[i][j];

        /* Reduce partial results stored in diffnorm by adding them.
           Leave final value in variable globaldiffnorm in all processes. */
	MPI_Iallreduce( &diffnorm, &globaldiffnorm, 1, MPI_DOUBLE, MPI_SUM,
		       MPI_COMM_WORLD, &r[0] ); /* Statement S19 */
	MPI_Wait( &r[0], &status ); /* Wait for the MPI_Iallreduce */
	globaldiffnorm = sqrt( globaldiffnorm );
	if (rank == 0) printf( "At iteration %d, diff is %e\n", itcnt, 
			       globaldiffnorm );
    } while (globaldiffnorm > 1.0e-2 && itcnt < 100);

    /* Collect into x the data segments distributed (xlocal) and print x. 
       Hint: each process sends its local segment (not the ghost areas)
             to process 0. */
    MPI_Igather( &xlocal[1], maxn * (maxn/size),MPI_DOUBLE,
		x, maxn * (maxn/size), MPI_DOUBLE, 
		0, MPI_COMM_WORLD, &r[0] ); /* Statement S20 */
    MPI_Wait( &r[0], &status ); /* Wait for the MPI_Igather */
    if (rank == 0) {
	printf( "Final solution is\n" );
	for (i=0; i<maxn; i++) {
	    for (j=0; j<maxn; j++) 
		printf( "%f ", x[i][j] );
	    printf( "\n" );
	}
    }

    MPI_Finalize( );
    return 0;
}
