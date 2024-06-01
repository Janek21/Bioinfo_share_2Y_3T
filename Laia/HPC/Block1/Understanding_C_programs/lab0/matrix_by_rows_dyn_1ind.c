/* C code which dynamically allocates memory for a 2 dimensional array (matrix).
   The 2D access to the matrix, normally done through indices i and j, is
   linearized into 1D. using a single index. Note how the row index i is
   multiplied by the dimension of the matrix (the number of columns since
   row-wise storage is used) to compute the offset from the beginning of the
   array to the beginning of row i. Then, the index of the column (j) is added
   to reach element in column j witin row i. */

#include <stdio.h>
#include <stdlib.h>

int main() {
    // Specify matrix dimensions
    const int rows = 3;
    const int cols = 3;

    // Dynamically allocate memory for the matrix
    int* matrix = (int*)malloc(rows * cols * sizeof(int));

    // Assign values to the matrix
    int count = 1;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix[i * cols + j] = count++;
        }
    }

    // Accessing matrix by rows
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            printf("%d ", matrix[i * cols + j]);
        }
        printf("\n");
    }

    // Deallocate dynamically allocated memory
    free(matrix);

    return 0;
}

