# Assuming matrix is a 2D list
matrix = [
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9]
]

# Accessing matrix by columns
for col in zip(*matrix):
    print(col)

