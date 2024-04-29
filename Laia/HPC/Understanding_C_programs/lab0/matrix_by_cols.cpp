/* Accessing a matrix by columns.
     Note that, in order to access the matrix by columns:
      - i indicates the row, since it is used as the first index when
          accessing the 2D array (matrix).
      - loop i is the inner loop.
*/
#include <iostream>
#include <vector>

int main() {
    // Assuming matrix is a vector of vectors
    std::vector<std::vector<int>> matrix = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    };

    // Accessing matrix by columns
    for (size_t j = 0; j < matrix[0].size(); ++j) {
        for (size_t i = 0; i < matrix.size(); ++i) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}

