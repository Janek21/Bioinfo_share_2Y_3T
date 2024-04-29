#include <iostream>
#include <vector>

int main() {
    // Assuming matrix is a vector of vectors
    std::vector<std::vector<int>> matrix = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    };

    // Accessing matrix by rows
    for (const auto &row : matrix) {
        for (int element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}

