#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#define _DEBUG 0 /* Define as 1 to get more information */

/* This code handles a number of processors (P) that does not divide the
   problem size (maxn) evenly. */

/* This example handles a maxn x maxn mesh, on P processors. */
/* P does not necessarily need to divide maxn evenly. */
#define maxn 13
#define P 5

/* Hint: Have a look at the auxiliary subroutine in the MPI matrix
         multiplication code in the slides. */

/* Adjust number of rows per process. */
int getRowCount(int rowsTotal, int mpiRank, int mpiSize) {
    /* Adjust slack of rows in case rowsTotal is not exactly divisible */
    return (rowsTotal / mpiSize) + mpiRank ; /* Statement S21 */
}


int main( argc, argv )
int argc;
char **argv;
{   
    int        rank, value, size, errcnt, toterr, i, j, itcnt;
    int        i_first, i_last;
    MPI_Status status;
    double     diffnorm, globaldiffnorm;
    double     xlocal[(maxn/P)+3][maxn];
    double     xnew[(maxn/P)+2][maxn];
    /* For the gatherv */
    double     x[maxn][maxn];
    int        lcnt;
    int        recvcnts[P];
    int        displs[P];
    int        nrows;

    MPI_Init( &argc, &argv );

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    /* xlocal[0][] is lower ghost area, xlocal[i_last+1][] is upper */
    /* Note that first and last processes have one less row of interior
       points */
    i_first = 1;
    /* Remove restriction to have maxn as multiple of P */
       // This was previously fixed: i_last  = maxn/size;
    nrows  = getRowCount(maxn, rank, size);; /* Statement S22 */
    i_last  = nrows;
    if (rank == 0)        i_first++;
    if (rank == size - 1) i_last--;
#if _DEBUG == 1
    printf(" Process %d: nrows=%d, i_first=%d i_last=%d\n",
             rank, nrows, i_first, i_last);
#endif

    /* Fill the data as specified */
    for (i=i_first; i<=i_last; i++) 
        for (j=0; j<maxn; j++) 
            xlocal[i][j] = rank;
    for (j=0; j<maxn; j++) {
        xlocal[i_first-1][j] = -1;
        xlocal[i_last+1][j] = -1;
    }

    itcnt = 0;
    do {
        /* Note the use of xlocal[i] for &xlocal[i][0] */

        /* Send last local row to next process (process with identifier rank+1)
           unless I'm the last one (process with identifier size-1). */
        if (rank < size - 1) 
            MPI_Send(xlocal[i_last],maxn,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD );
        /* Receive from previous unless I'm the first process (rank==0)
           and store row in initial row (ghost area) */
        if (rank > 0)
            MPI_Recv(xlocal[0],maxn,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD, &status );
        /* Send 1st local row to previous process unless I'm the 1st process */
        if (rank > 0) 
            MPI_Send(xlocal[1],maxn,MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD );
        /* Receive from next process and store row in last row in xlocal
           (ghost area) unless I'm the last process */
        if (rank < size - 1) 
            MPI_Recv(xlocal[i_last+1],maxn,MPI_DOUBLE,rank+1,1,MPI_COMM_WORLD, &status );
        
        /* Compute new values (but not on boundary) */
        itcnt ++;
        diffnorm = 0.0;
        for (i=i_first; i<=i_last; i++) {
            //printf(" Process %d: Loop Compute new values i=%d\n", rank, i);
            for (j=1; j<maxn-1; j++) {
                xnew[i][j] = (xlocal[i][j+1] + xlocal[i][j-1] +
                              xlocal[i+1][j] + xlocal[i-1][j]) / 4.0;
                diffnorm += (xnew[i][j] - xlocal[i][j]) * 
                            (xnew[i][j] - xlocal[i][j]);
            }
        }
        /* Only transfer the interior points */
        //for (i=i_first; i<=i_last; i++) 
        for (i=i_first; i<=i_last; i++) {
            //printf(" Process %d: Loop transfer segment i=%d\n", rank, i);
            for (j=1; j<maxn-1; j++) 
                xlocal[i][j] = xnew[i][j];
        }

        MPI_Allreduce( &diffnorm,&globaldiffnorm,1,MPI_DOUBLE,MPI_SUM,
                       MPI_COMM_WORLD );
        globaldiffnorm = sqrt( globaldiffnorm );
//      if (rank == 0 ) printf( "At iteration %d, diff is %e\n", itcnt, 
//                             globaldiffnorm );
    } while (globaldiffnorm > 1.0e-2 && itcnt < 200);
    if (rank == 0 ) printf( "At iteration %d, diff is %e\n", itcnt, 
                               globaldiffnorm );
#if _DEBUG == 2
    printf(" Process %d: Last row (xlocal[%d]) in local segment\n", rank, i_last);
            for (j=0; j<maxn; j++) 
                printf( "%f ", xlocal[i_last][j] );
            printf( "\n" );
#endif

    /* Collect the data into x and print it */
    /* Use a preliminary gather to get the recv counts */
    /* Next statement was
         lcnt = maxn * (maxn / size);
       Adjusted to account for any number of rows in different processes */
    lcnt = maxn * nrows; /* Total number of values in local segment */
    /* Inform the master process of the number of values in local segment */ 
    MPI_Gather( &lcnt, 0 , MPI_INT , recvcnts, 1, MPI_INT, 0 , MPI_COMM_WORLD ); /* Statement S23 */
    /* Form the displacements using the recv counts just obtained */
    displs[0] = 0;
    for (i=1; i<size; i++) 
        displs[i] = displs[i-1] + recvcnts[i-1];

    /* Now gather with proper amount of values (lcnt) from each process */
    MPI_Gatherv( xlocal[1], 1 , MPI_DOUBLE, x, 0 , displs, MPI_DOUBLE, 0, MPI_COMM_WORLD ); /* Statement S24 */
    if (rank == 0) {
#if _DEBUG == 1
        for (i=0; i<size; i++) 
            printf("recvcnts[%d]=%d\tdispls[%d]=%d\n",
                    i, recvcnts[i], i, displs[i]);
#endif
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

