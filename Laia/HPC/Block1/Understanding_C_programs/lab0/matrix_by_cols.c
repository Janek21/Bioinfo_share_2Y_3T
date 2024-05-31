/* Accessing a matrix by columns.
     Note that, in order to access the matrix by columns:
      - i indicates the row, since it is used as the first index when
          accessing the 2D array (matrix).
      - loop i is the inner loop.
*/
#include <stdio.h>

int main() {
    // Assuming matrix is a 2D array stored by rows.
    int matrix[3][3] = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    };

    // Accessing matrix by columns
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }

    return 0;
}

