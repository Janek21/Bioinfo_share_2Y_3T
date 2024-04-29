/* C++ code: dynamically allocates memory for a 2 dimensional array (matrix).
   The 2D access to the matrix, normally done through indices i and j, is
   linearized into 1D. using a single index. Note how the row index i is
   multiplied by the dimension of the matrix (the number of columns since
   row-wise storage is used) to compute the offset from the beginning of the
   array to the beginning of row i. Then, the index of the column (j) is added
   to reach element in column j witin row i. */

#include <iostream>

int main() {
    // Specify matrix dimensions
    const int rows = 3;
    const int cols = 3;

    // Dynamically allocate memory for the matrix
    int* matrix = new int[rows * cols];

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
            std::cout << matrix[i * cols + j] << " ";
        }
        std::cout << std::endl;
    }

    // Deallocate dynamically allocated memory
    delete[] matrix;

    return 0;
}

