#include <stdio.h>
#include <stdlib.h>

int main() {
    // Specify matrix dimensions
    const int rows = 3;
    const int cols = 3;

    // Dynamically allocate memory for the matrix
    int** matrix = (int**)malloc(rows * sizeof(int*));
    for (int i = 0; i < rows; ++i) {
        matrix[i] = (int*)malloc(cols * sizeof(int));
    }

    // Assign values to the matrix
    int count = 1;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix[i][j] = count++;
        }
    }

    // Accessing matrix by rows
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }

    // Deallocate dynamically allocated memory
    for (int i = 0; i < rows; ++i) {
        free(matrix[i]);
    }
    free(matrix);

    return 0;
}

