
/* Column-wise storage in C

   If you  want to  store the  matrix by columns  and keep  the column-wise
   access,  you need  to adjust  the  memory allocation  and matrix  access
   accordingly. Note tha the matrix access  is now linearized into 1D using
   the expression  "j * rows +  i", instead of "i  * cols + j"  as was done
   assuming row-wise storage. 
*/

#include <stdio.h>
#include <stdlib.h>

int main() {
    // Specify matrix dimensions
    const int rows = 3;
    const int cols = 3;

    // Dynamically allocate memory for the matrix stored by columns
    int* matrix = (int*)malloc(rows * cols * sizeof(int));

    // Assign values to the matrix
    int count = 1;
    for (int j = 0; j < cols; ++j) {
        for (int i = 0; i < rows; ++i) {
            matrix[j * rows + i] = count++;
        }
    }

    // Accessing matrix by columns
    for (int j = 0; j < cols; ++j) {
        printf("Column %d: (", j);
        for (int i = 0; i < rows; ++i) {
            printf("%d ", matrix[j * rows + i]);
        }
        printf(")'\n");
    }
        
    printf("Note: (. . .)' stands for the transposition of the vector.\n");

    // Deallocate dynamically allocated memory
    free(matrix);

    return 0;
}
